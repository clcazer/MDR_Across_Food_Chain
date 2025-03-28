
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

source(here("Scripts", "function_scripts", "data_wrangling.R"))
source(here("Scripts", "function_scripts", "rule_mining_and_selection.R"))
source(here("Scripts", "function_scripts", "rule_analysis.R"))
source(here("Scripts", "function_scripts", "network_graphs.R"))






#get the data ready
cecal_2017Onward_phenotype_df <- read.csv("dataset_specific_outputs/cecal_2017Onward/NEW_cecal_2017Onward_wide_resStatus_phenotype.csv")
cecal_2017Onward_Class_level_phenotype_df <- read.csv("dataset_specific_outputs/cecal_2017Onward/NEW_cecal_2017Onward_wide_class_level_phenotype.csv")



cecal_2017Onward_genotype_df <- read.csv("dataset_specific_outputs/cecal_2017Onward/cecal_2017Onward_wide_corrected_genotype.csv")
cecal_2017Onward_Family_genotype_df <- read.csv("dataset_specific_outputs/cecal_2017Onward/cecal_2017Onward_wide_family_level_genotype.csv")
cecal_2017Onward_Class_level_genotype_df <- read.csv("dataset_specific_outputs/cecal_2017Onward/cecal_2017Onward_wide_class_level_genotype.csv")



run_rule_mining_and_selection <- function(df, class_df = NULL, data_source, resistance_indicator) {

print("RULE MINING AND SELECTION STARTED")

#     #get histogram of top pc1, pc2, pc3, and pc4 loadings
#      select_quality_measures(plothist = TRUE, resistance_indicator = resistance_indicator, 
#      df = df, target = "rules", data_source = data_source)
# print("1")
#      select_measure_cutoffs(selected_measures = c("cosine", "jaccard", "kulczynski", "support"), 
#      df = df, resistance_indicator = resistance_indicator, target = "rules", data_source = data_source)

# print("2")
#      rarefaction_by_sample_size(df = df, resistance_indicator = resistance_indicator, target = "rules", data_source = data_source)

# print("3")

#     rarefaction_rule_v_class(df = df, resistance_indicator = resistance_indicator, target = "rules", data_source = data_source)

# print("4")

#     unique_classes_represented(df = df, resistance_indicator = resistance_indicator, target = "rules", data_source = data_source)

# print("5")

#     rules_v_cutoffs(df = df, resistance_indicator = resistance_indicator, target = "rules", measures = c("cosine", "jaccard", "kulczynski", "support"), low = 0, high = 1, data_source = data_source)

# print("6")

    prevalence_descriptives(df = df, data_source = data_source, resistance_indicator = resistance_indicator)

print("7")
    prevalence_descriptives(df = class_df, data_source = data_source, resistance_indicator = resistance_indicator, class_level = TRUE)
  

print("RULE MINING AND SELECTION DONE")
}


# RUN THE FOLLOWING TWO LINES TO GET THE RULE MINING AND SELECTION DONE
run_rule_mining_and_selection(df = cecal_2017Onward_phenotype_df, class_df = cecal_2017Onward_Class_level_phenotype_df, data_source = "cecal_2017Onward", resistance_indicator = "phenotype")
run_rule_mining_and_selection(df = cecal_2017Onward_genotype_df, class_df = cecal_2017Onward_Class_level_genotype_df, data_source = "cecal_2017Onward", resistance_indicator = "genotype")


#plot all vs best rules
plot_all_vs_best_rules(df = cecal_2017Onward_phenotype_df, resistance_indicator = "phenotype", target = "rules", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward")
plot_all_vs_best_rules(df = cecal_2017Onward_Family_genotype_df, resistance_indicator = "genotype", target = "rules", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward")

plot_all_vs_best_rules(df = cecal_2017Onward_Class_level_phenotype_df, resistance_indicator = "phenotype", target = "rules", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward", class_level = TRUE)
plot_all_vs_best_rules(df = cecal_2017Onward_Class_level_genotype_df, resistance_indicator = "genotype", target = "rules", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward", class_level = TRUE)





#NETWORK GRAPHS WERE NOT USED IN PAPER
#create Network Graphs

# #phenotype (not class level)
#  graph_rules(df = cecal_2017Onward_phenotype_df, target = "rules", cut_off = c(0.5, 0, 0.5, 0), resistance_indicator = "phenotype", measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward", rules_selected = "best", agg = FALSE)
# #genotype (not class level)
#  graph_rules(df = cecal_2017Onward_Family_genotype_df, target = "rules", cut_off = c(0.5, 0, 0.5, 0), resistance_indicator = "genotype", measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward", rules_selected = "best", agg = FALSE)

# #phenotype (class level)
#  graph_rules(df = cecal_2017Onward_Class_level_phenotype_df, target = "rules", cut_off = c(0.5, 0, 0.5, 0), resistance_indicator = "phenotype", measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward", rules_selected = "best", agg = TRUE)
# #genotype (class level)
#  graph_rules(df = cecal_2017Onward_Class_level_genotype_df, target = "rules", cut_off = c(0.5, 0, 0.5, 0), resistance_indicator = "genotype", measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward", rules_selected = "best", agg = TRUE)





print("DONE RUNNING CECAL 2017 ONWARD ANALYSIS")