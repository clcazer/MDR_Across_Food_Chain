
library(stringr)
library(rlist)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(readxl)
library(cowplot)
library(arules)
library(plyr)
library(dplyr)

library(here)

source(here("function_scripts", "data_wrangling.R"))
source(here("function_scripts", "rule_mining_and_selection.R"))
source(here("function_scripts", "rule_analysis.R"))
source(here("function_scripts", "network_graphs.R"))




#load in the cecal data
cecal_phenotype_df <- read.csv("cecal/cecal_wide_resStatus_phenotype.csv")
cecal_genotype_df <- read.csv("cecal/cecal_wide_family_level_genotype.csv")


cecal_Class_level_phenotype_df <- read.csv("cecal/cecal_wide_class_level_phenotype.csv")
cecal_Class_level_genotype_df <- read.csv("cecal/cecal_wide_class_level_genotype.csv")


#add a "cecal" column and fill with 1's
cecal_phenotype_df$cecal <- 1
cecal_genotype_df$cecal <- 1
cecal_Class_level_phenotype_df$cecal <- 1
cecal_Class_level_genotype_df$cecal <- 1
#add a Retail_Meats column and fill with 0's
cecal_phenotype_df$Retail_Meats <- 0
cecal_genotype_df$Retail_Meats <- 0
cecal_Class_level_phenotype_df$Retail_Meats <- 0
cecal_Class_level_genotype_df$Retail_Meats <- 0



#load in the Retail_Meats data
Retail_Meats_phenotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_resStatus_phenotype.csv")
Retail_Meats_genotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_family_level_genotype.csv")    

Retail_Meats_Class_level_phenotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_phenotype.csv")
Retail_Meats_Class_level_genotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_genotype.csv")

#add a "cecal" column and fill with 0's
Retail_Meats_phenotype_df$cecal <- 0
Retail_Meats_genotype_df$cecal <- 0
Retail_Meats_Class_level_phenotype_df$cecal <- 0
Retail_Meats_Class_level_genotype_df$cecal <- 0
#add a "Retail_Meats" column and fill with 1's
Retail_Meats_phenotype_df$Retail_Meats <- 1
Retail_Meats_genotype_df$Retail_Meats <- 1
Retail_Meats_Class_level_phenotype_df$Retail_Meats <- 1
Retail_Meats_Class_level_genotype_df$Retail_Meats <- 1

#for both cecal and retail meats dfs, only include years that are present in both cecal and retail meats
cecal_phenotype_df <- cecal_phenotype_df[cecal_phenotype_df$Year %in% Retail_Meats_phenotype_df$Year, ]
cecal_genotype_df <- cecal_genotype_df[cecal_genotype_df$Year %in% Retail_Meats_genotype_df$Year, ]
Retail_Meats_phenotype_df <- Retail_Meats_phenotype_df[Retail_Meats_phenotype_df$Year %in% cecal_phenotype_df$Year, ]
Retail_Meats_genotype_df <- Retail_Meats_genotype_df[Retail_Meats_genotype_df$Year %in% cecal_genotype_df$Year, ]
#do this for the class level dfs as well
cecal_Class_level_phenotype_df <- cecal_Class_level_phenotype_df[cecal_Class_level_phenotype_df$Year %in% Retail_Meats_Class_level_phenotype_df$Year, ]
cecal_Class_level_genotype_df <- cecal_Class_level_genotype_df[cecal_Class_level_genotype_df$Year %in% Retail_Meats_Class_level_genotype_df$Year, ]
Retail_Meats_Class_level_phenotype_df <- Retail_Meats_Class_level_phenotype_df[Retail_Meats_Class_level_phenotype_df$Year %in% cecal_Class_level_phenotype_df$Year, ]
Retail_Meats_Class_level_genotype_df <- Retail_Meats_Class_level_genotype_df[Retail_Meats_Class_level_genotype_df$Year %in% cecal_Class_level_genotype_df$Year, ]

print("PREPROCESSED THE DATA")


# Function to add missing columns with zeros
add_missing_columns <- function(df1, df2) {
  missing_cols <- setdiff(names(df2), names(df1))
  if (length(missing_cols) > 0) {
    df1[missing_cols] <- 0
  }
  return(df1)
}

# Add missing columns to each dataframe before combining
Retail_Meats_phenotype_df <- add_missing_columns(Retail_Meats_phenotype_df, cecal_phenotype_df)
cecal_phenotype_df <- add_missing_columns(cecal_phenotype_df, Retail_Meats_phenotype_df)

Retail_Meats_genotype_df <- add_missing_columns(Retail_Meats_genotype_df, cecal_genotype_df)
cecal_genotype_df <- add_missing_columns(cecal_genotype_df, Retail_Meats_genotype_df)

Retail_Meats_Class_level_phenotype_df <- add_missing_columns(Retail_Meats_Class_level_phenotype_df, cecal_Class_level_phenotype_df)
cecal_Class_level_phenotype_df <- add_missing_columns(cecal_Class_level_phenotype_df, Retail_Meats_Class_level_phenotype_df)

Retail_Meats_Class_level_genotype_df <- add_missing_columns(Retail_Meats_Class_level_genotype_df, cecal_Class_level_genotype_df)
cecal_Class_level_genotype_df <- add_missing_columns(cecal_Class_level_genotype_df, Retail_Meats_Class_level_genotype_df)

print("ADDED MISSING COLUMNS")
#combine Retail_meats and cecal phenotype dataframes
combined_phenotype_df <- rbind(Retail_Meats_phenotype_df, cecal_phenotype_df)
#combine Retail_meats and cecal genotype dataframes
combined_genotype_df <- rbind(Retail_Meats_genotype_df, cecal_genotype_df)
#combine class level phenotype dataframes
combined_Class_level_phenotype_df <- rbind(Retail_Meats_Class_level_phenotype_df, cecal_Class_level_phenotype_df)
#combine class level genotype dataframes
combined_Class_level_genotype_df <- rbind(Retail_Meats_Class_level_genotype_df, cecal_Class_level_genotype_df)









run_rule_mining_and_selection <- function(df, class_df = NULL, data_source, resistance_indicator) {

print("RULE MINING AND SELECTION STARTED")

    #get histogram of top pc1, pc2, pc3, and pc4 loadings
     select_quality_measures(plothist = TRUE, resistance_indicator = resistance_indicator, 
     df = df, target = "rules", data_source = data_source)
print("1")
     select_measure_cutoffs(selected_measures = c("cosine", "jaccard", "kulczynski", "support"), 
     df = df, resistance_indicator = resistance_indicator, target = "rules", data_source = data_source)

print("2")
     rarefaction_by_sample_size(df = df, resistance_indicator = resistance_indicator, target = "rules", data_source = data_source)

print("3")

    rarefaction_rule_v_class(df = df, resistance_indicator = resistance_indicator, target = "rules", data_source = data_source)

print("4")

    unique_classes_represented(df = df, resistance_indicator = resistance_indicator, target = "rules", data_source = data_source)

print("5")

    rules_v_cutoffs(df = df, resistance_indicator = resistance_indicator, target = "rules", measures = c("cosine", "jaccard", "kulczynski", "support"), low = 0, high = 1, data_source = data_source)

print("6")

    prevalence_descriptives(df = df, data_source = data_source, resistance_indicator = resistance_indicator)

print("7")
    prevalence_descriptives(df = class_df, data_source = data_source, resistance_indicator = resistance_indicator, class_level = TRUE)
  

print("RULE MINING AND SELECTION DONE")
}


# RUN THE FOLLOWING TWO LINES TO GET THE RULE MINING AND SELECTION DONE
#run_rule_mining_and_selection(df = combined_phenotype_df, class_df = combined_Class_level_phenotype_df, data_source = "combined_source_encodings", resistance_indicator = "phenotype")
#run_rule_mining_and_selection(df = combined_genotype_df, class_df = combined_Class_level_genotype_df, data_source = "combined_source_encodings", resistance_indicator = "genotype")







#create Network Graphs


#phenotype (not class level)
 graph_rules(df = combined_phenotype_df, target = "rules", cut_off = c(0.5, 0, 0.5, 0), resistance_indicator = "phenotype", measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "combined_source_encodings", rules_selected = "all", agg = FALSE)
#genotype (not class level)
 #graph_rules(df = combined_genotype_df, target = "rules", cut_off = c(0.5, 0, 0.5, 0), resistance_indicator = "genotype", measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "combined_source_encodings", rules_selected = "best", agg = FALSE)

#phenotype (class level)
 graph_rules(df = combined_Class_level_phenotype_df, target = "rules", cut_off = c(0.5, 0, 0.5, 0), resistance_indicator = "phenotype", measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "combined_source_encodings", rules_selected = "all", agg = TRUE)
#genotype (class level)
 graph_rules(df = combined_Class_level_genotype_df, target = "rules", cut_off = c(0.5, 0, 0.5, 0), resistance_indicator = "genotype", measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "combined_source_encodings", rules_selected = "all", agg = TRUE)



#plot all vs best rules
plot_all_vs_best_rules(df = combined_phenotype_df, resistance_indicator = "phenotype", target = "rules", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "combined_source_encodings")
plot_all_vs_best_rules(df = combined_genotype_df, resistance_indicator = "genotype", target = "rules", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "combined_source_encodings")

plot_all_vs_best_rules(df = combined_Class_level_phenotype_df, resistance_indicator = "phenotype", target = "rules", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "combined_source_encodings", class_level = TRUE)
plot_all_vs_best_rules(df = combined_Class_level_genotype_df, resistance_indicator = "genotype", target = "rules", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "combined_source_encodings", class_level = TRUE)


