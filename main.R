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
source("data_wrangling.R")
source("rule_mining_and_selection.R")
source("rule_analysis.R")






#function to run all the rule mining and selection for a given (individual) dataset (data_source) and resistance indicator (phenotype or genotype)

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

cecal_phenotype_df <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_resStatus_phenotype.csv")
cecal_Class_phenotype_df <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_class_level_phenotype.csv")
cecal_Family_genotype_df <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_family_level_genotype.csv")
cecal_Class_genotype_df <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_class_level_genotype.csv")

Retail_Meats_phenotype_df <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_resStatus_phenotype.csv")
Retail_Meats_Class_phenotype_df <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_class_level_phenotype.csv")
Retail_Meats_Family_genotype_df <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_family_level_genotype.csv")
Retail_Meats_Class_genotype_df <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_class_level_genotype.csv")

#   run_rule_mining_and_selection(df = Retail_Meats_phenotype_df, class_df = Retail_Meats_Class_phenotype_df, data_source = "Retail_Meats_2017Onward", resistance_indicator = "phenotype")
#   run_rule_mining_and_selection(df = Retail_Meats_Family_genotype_df, class_df = Retail_Meats_Class_genotype_df, data_source = "Retail_Meats_2017Onward", resistance_indicator = "genotype")
#   run_rule_mining_and_selection(df = cecal_phenotype_df, class_df = cecal_Class_phenotype_df, data_source = "cecal_2017Onward", resistance_indicator = "phenotype")
#   run_rule_mining_and_selection(df = cecal_Family_genotype_df, class_df = cecal_Class_genotype_df, data_source = "cecal_2017Onward", resistance_indicator = "genotype")









run_rule_analysis_within_indicators <- function(df, class_df = NULL, data_source, resistance_indicator, rules_selected, class_level = FALSE) {
        print("RULE ANALYSIS WITHIN INDICATORS STARTED")
        
        #plot all vs best rules
        plot_all_vs_best_rules(df = df, resistance_indicator = resistance_indicator, target = "rules", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = data_source)
        print("1")
        #do rule overlap and crs
        rule_overlap(df = df, resistance_indicator = resistance_indicator, target = "rules", rules_selected = rules_selected, cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = data_source)
        print("2")
        cumulative_rule_stability(df = df, resistance_indicator = resistance_indicator, target = "rules", rules_selected = rules_selected, cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = data_source)
        print("3")
        rule_overlap(df = class_df, class_level = TRUE, resistance_indicator = resistance_indicator, target = "rules", rules_selected = rules_selected, cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = data_source)
        print("4")
        cumulative_rule_stability(df = class_df, class_level = TRUE, resistance_indicator = resistance_indicator, target = "rules", rules_selected = rules_selected, cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = data_source)
        print("5")

        #do proportion captured and cumulative proportion captured
        proportion_captured(df = df, resistance_indicator = resistance_indicator, target = "rules", rules_selected = rules_selected, cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = data_source)
        print("6")
        cumulative_proportion_captured(df = df, resistance_indicator = resistance_indicator, target = "rules", rules_selected = rules_selected, cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = data_source)
        print("7")
        proportion_captured(df = class_df, class_level = TRUE, resistance_indicator = resistance_indicator, target = "rules", rules_selected = rules_selected, cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = data_source)
        print("8")
        cumulative_proportion_captured(df = class_df, class_level = TRUE, resistance_indicator = resistance_indicator, target = "rules", rules_selected = rules_selected, cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = data_source)


print("RULE ANALYSIS WITHIN INDICATORS DONE")

}


#   run_rule_analysis_within_indicators(df = Retail_Meats_phenotype_df, class_df = Retail_Meats_Class_phenotype_df, data_source = "Retail_Meats_2017Onward", resistance_indicator = "phenotype", rules_selected = "best")
#   run_rule_analysis_within_indicators(df = Retail_Meats_Family_genotype_df, class_df = Retail_Meats_Class_genotype_df, data_source = "Retail_Meats_2017Onward", resistance_indicator = "genotype", rules_selected = "best")
#   run_rule_analysis_within_indicators(df = cecal_phenotype_df, class_df = cecal_Class_phenotype_df, data_source = "cecal_2017Onward", resistance_indicator = "phenotype", rules_selected = "best")
#   run_rule_analysis_within_indicators(df = cecal_Family_genotype_df, class_df = cecal_Class_genotype_df, data_source = "cecal_2017Onward", resistance_indicator = "genotype", rules_selected = "best")



run_rule_analysis_between_indicators <- function(genotype_df, phenotype_df, target, rules_selected, cut_off, measures_used, data_source, reference_indicator, comparison_indicator) {


        gene_vs_phenotype_cumulative_percent_captured(genotype_df, phenotype_df, target, rules_selected, cut_off, measures_used, data_source, reference_indicator, comparison_indicator)
        
}

# Include only rows where the ID matches in both dataframes
Retail_Meats_Class_phenotype_df <- Retail_Meats_Class_phenotype_df[Retail_Meats_Class_phenotype_df$ID %in% Retail_Meats_Class_genotype_df$ID, ]
Retail_Meats_Class_genotype_df <- Retail_Meats_Class_genotype_df[Retail_Meats_Class_genotype_df$ID %in% Retail_Meats_Class_phenotype_df$ID, ]


cecal_Class_phenotype_df <- cecal_Class_phenotype_df[cecal_Class_phenotype_df$ID %in% cecal_Class_genotype_df$ID, ]
cecal_Class_genotype_df <- cecal_Class_genotype_df[cecal_Class_genotype_df$ID %in% cecal_Class_phenotype_df$ID, ]

run_rule_analysis_between_indicators(genotype_df = Retail_Meats_Class_genotype_df, phenotype_df = Retail_Meats_Class_phenotype_df, target = "rules", rules_selected = "all", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "Retail_Meats_2017Onward", reference_indicator = "Phenotype", comparison_indicator = "Genotype")
run_rule_analysis_between_indicators(genotype_df = Retail_Meats_Class_genotype_df, phenotype_df = Retail_Meats_Class_phenotype_df, target = "rules", rules_selected = "all", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "Retail_Meats_2017Onward", reference_indicator = "Genotype", comparison_indicator = "Phenotype")



 run_rule_analysis_between_indicators(genotype_df = cecal_Class_genotype_df, phenotype_df = cecal_Class_phenotype_df, target = "rules", rules_selected = "all", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward", reference_indicator = "Phenotype", comparison_indicator = "Genotype")
 run_rule_analysis_between_indicators(genotype_df = cecal_Class_genotype_df, phenotype_df = cecal_Class_phenotype_df, target = "rules", rules_selected = "all", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "cecal_2017Onward", reference_indicator = "Genotype", comparison_indicator = "Phenotype")












