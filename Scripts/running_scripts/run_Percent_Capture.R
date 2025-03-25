
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

source(here("Scripts","function_scripts", "data_wrangling.R"))
source(here("Scripts","function_scripts", "rule_mining_and_selection.R"))
source(here("Scripts","function_scripts", "rule_analysis.R"))



#####################################################################################################################################################
#####################################################################################################################################################
################################### DEFINE A FUNCTION TO RUN ANALYSIS FOR EACH DATA PARTITION #######################################################
#####################################################################################################################################################
#####################################################################################################################################################

run_all_percent_captured_analyses <- function(data_partition,

                                            cecal_phenotype_df, 
                                            cecal_Class_phenotype_df,

                                            cecal_Family_genotype_df ,
                                            cecal_Class_genotype_df,

                                            Retail_Meats_phenotype_df,
                                            Retail_Meats_Class_phenotype_df,

                                            Retail_Meats_Family_genotype_df,
                                            Retail_Meats_Class_genotype_df,




                                            cecal_2017_phenotype_df, 
                                            cecal_2017_Class_phenotype_df,

                                            cecal_2017_Family_genotype_df ,
                                            cecal_2017_Class_genotype_df,

                                            Retail_2017_Meats_phenotype_df,
                                            Retail_2017_Meats_Class_phenotype_df,

                                            Retail_2017_Meats_Family_genotype_df,
                                            Retail_2017_Meats_Class_genotype_df,



                                            NAHLN_phenotype_df, 
                                            NAHLN_Class_phenotype_df,

                                            NAHLN_Family_genotype_df ,
                                            NAHLN_Class_genotype_df


                                                ) {

    #####################################
#RUN THE WITHIN INDICATORS COMPARISON
#####################################

print("STARTED WITHIN INDICATORS PLOT")

phenotype_datasets <- list(
    "Cecal Phenotype" = cecal_phenotype_df,
    "Retail Meats Phenotype" = Retail_Meats_phenotype_df,
    "NAHLN Phenotype" = NAHLN_phenotype_df
)

class_phenotype_datasets <- list(
    "Cecal Class Phenotype" = cecal_Class_phenotype_df,
    "Retail Meats Class Phenotype" = Retail_Meats_Class_phenotype_df,
    "NAHLN Class Phenotype" = NAHLN_Class_phenotype_df
)

genotype_datasets <- list(
    "Cecal Genotype" = cecal_Family_genotype_df,
    "Retail Meats Genotype" = Retail_Meats_Family_genotype_df,
    "NAHLN Genotype" = NAHLN_Family_genotype_df
)

class_genotype_datasets <- list(
    "Cecal Class Genotype" = cecal_Class_genotype_df,
    "Retail Meats Class Genotype" = Retail_Meats_Class_genotype_df,
    "NAHLN Class Genotype" = NAHLN_Class_genotype_df
)

 nonclass_plot <- cumulative_proportion_captured_multi_v2(
    datasets_geno = genotype_datasets,
    datasets_pheno = phenotype_datasets,
    class_level = FALSE,
    target = "rules",
    rules_selected = "best",
    cut_off = c(0.5, 0, 0.5, 0),
    measures_used = c("cosine", "jaccard", "kulczynski", "support"),
    data_partition = data_partition
)

 class_plot <- cumulative_proportion_captured_multi_v2(
    datasets_geno = class_genotype_datasets,
    datasets_pheno = class_phenotype_datasets,
    class_level = TRUE,
    target = "rules",
    rules_selected = "best",
    cut_off = c(0.5, 0, 0.5, 0),
    measures_used = c("cosine", "jaccard", "kulczynski", "support"),
    data_partition = data_partition
)






format_and_combine_plots <- function(plot_list, title_text, labels, ncol = 2, nrow =1, title = FALSE, legend_size = 12) {
    # Format individual plots
    formatted_plots <- lapply(seq_along(plot_list), function(i) {
        p <- plot_list[[i]]
        p <- p + 
            theme(
                plot.title = element_blank(),  # Remove individual plot titles
                axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
                axis.title.x = element_text(margin = margin(t = 10)),
                axis.title.y = element_text(margin = margin(r = 10)),
                legend.text = element_text(size = legend_size),
                legend.key.width = unit(2, "cm")
            ) +
            scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1))
        
        # Remove y-axis title for right column plots
        if (i %% 2 == 0) {
            p <- p + theme(axis.title.y = element_blank())
        }
        
        return(p)
    })
    
    # Create panel of plots without title
    combined_plot <- plot_grid(
        plotlist = formatted_plots,
        ncol = ncol,
        nrow = nrow,
        labels = labels,
        label_size = 25,
        align = "hv",
        axis = "bt",
        rel_widths = c(1.05, 1)
    )

    return(combined_plot)
}

within_indicators_plot <- format_and_combine_plots(
    plot_list = list(nonclass_plot, class_plot),
    title_text = "",  # Remove title
    labels = c("A", "B")
)

print("DONE WITH WITHIN INDICATORS PLOT")
# Save the final plot with increased dimensions
# ggsave(plot = within_indicators_plot, filename = str_glue("combined_figures_outputs/percent_captured_analysis/within_indicators/rule_set_comparison_within_indicators.png"), 
#        width = 18, height = 10, dpi = 300)


#####################################
#RUN THE BETWEEN INDICATORS COMPARISON
#####################################
print("STARTED BETWEEN INDICATORS PLOT")

#since we're comparing between indicators, data has to be at the class level

#Include only rows where the ID matches in both dataframes for each data source
IDMathced_Retail_Meats_Class_phenotype_df <- Retail_Meats_Class_phenotype_df[Retail_Meats_Class_phenotype_df$ID %in% Retail_Meats_Class_genotype_df$ID, ]
IDMathced_Retail_Meats_Class_genotype_df <- Retail_Meats_Class_genotype_df[Retail_Meats_Class_genotype_df$ID %in% IDMathced_Retail_Meats_Class_phenotype_df$ID, ]

IDMathced_cecal_Class_phenotype_df <- cecal_Class_phenotype_df[cecal_Class_phenotype_df$ID %in% cecal_Class_genotype_df$ID, ]
IDMathced_cecal_Class_genotype_df <- cecal_Class_genotype_df[cecal_Class_genotype_df$ID %in% IDMathced_cecal_Class_phenotype_df$ID, ]

# Repeat the process for the 2017 data frames
IDMathced_Retail_2017_Meats_Class_phenotype_df <- Retail_2017_Meats_Class_phenotype_df[Retail_2017_Meats_Class_phenotype_df$ID %in% Retail_2017_Meats_Class_genotype_df$ID, ]
IDMathced_Retail_2017_Meats_Class_genotype_df <- Retail_2017_Meats_Class_genotype_df[Retail_2017_Meats_Class_genotype_df$ID %in% IDMathced_Retail_2017_Meats_Class_phenotype_df$ID, ]

IDMathced_cecal_2017_Class_phenotype_df <- cecal_2017_Class_phenotype_df[cecal_2017_Class_phenotype_df$ID %in% cecal_2017_Class_genotype_df$ID, ]
IDMathced_cecal_2017_Class_genotype_df <- cecal_2017_Class_genotype_df[cecal_2017_Class_genotype_df$ID %in% IDMathced_cecal_2017_Class_phenotype_df$ID, ]

# Get common columns between genotype and phenotype dataframes for each source
retail_common_cols <- intersect(colnames(IDMathced_Retail_Meats_Class_genotype_df), 
                              colnames(IDMathced_Retail_Meats_Class_phenotype_df))
cecal_common_cols <- intersect(colnames(IDMathced_cecal_Class_genotype_df),
                             colnames(IDMathced_cecal_Class_phenotype_df))

# Keep only common columns for retail dataframes
IDMathced_Retail_Meats_Class_genotype_df <- IDMathced_Retail_Meats_Class_genotype_df[, retail_common_cols]
IDMathced_Retail_Meats_Class_phenotype_df <- IDMathced_Retail_Meats_Class_phenotype_df[, retail_common_cols]

# Keep only common columns for cecal dataframes  
IDMathced_cecal_Class_genotype_df <- IDMathced_cecal_Class_genotype_df[, cecal_common_cols]
IDMathced_cecal_Class_phenotype_df <- IDMathced_cecal_Class_phenotype_df[, cecal_common_cols]

# Repeat for 2017 data
retail_2017_common_cols <- intersect(colnames(IDMathced_Retail_2017_Meats_Class_genotype_df),
                                   colnames(IDMathced_Retail_2017_Meats_Class_phenotype_df))
cecal_2017_common_cols <- intersect(colnames(IDMathced_cecal_2017_Class_genotype_df),
                                  colnames(IDMathced_cecal_2017_Class_phenotype_df))

IDMathced_Retail_2017_Meats_Class_genotype_df <- IDMathced_Retail_2017_Meats_Class_genotype_df[, retail_2017_common_cols]
IDMathced_Retail_2017_Meats_Class_phenotype_df <- IDMathced_Retail_2017_Meats_Class_phenotype_df[, retail_2017_common_cols]

IDMathced_cecal_2017_Class_genotype_df <- IDMathced_cecal_2017_Class_genotype_df[, cecal_2017_common_cols]
IDMathced_cecal_2017_Class_phenotype_df <- IDMathced_cecal_2017_Class_phenotype_df[, cecal_2017_common_cols]





# Create the dataset lists
genotype_datasets <- list(
    "Cecal Genotype" = IDMathced_cecal_Class_genotype_df,
    "Retail Meats Genotype" = IDMathced_Retail_Meats_Class_genotype_df,
    "NAHLN Genotype" = NAHLN_Class_genotype_df
)

phenotype_datasets <- list(
    "Cecal Phenotype" = IDMathced_cecal_Class_phenotype_df,
    "Retail Meats Phenotype" = IDMathced_Retail_Meats_Class_phenotype_df,
    "NAHLN Phenotype" = NAHLN_Class_phenotype_df
)

genotype2017_datasets <- list(
    "Cecal Genotype" = IDMathced_cecal_2017_Class_genotype_df,
    "Retail Meats Genotype" = IDMathced_Retail_2017_Meats_Class_genotype_df
)

phenotype2017_datasets <- list(
    "Cecal Phenotype" = IDMathced_cecal_2017_Class_phenotype_df,
    "Retail Meats Phenotype" = IDMathced_Retail_2017_Meats_Class_phenotype_df
)




indicator_comparison_plot <- gene_vs_phenotype_cumulative_percent_captured_multi_v2(
     datasets_geno = genotype_datasets,
     datasets_pheno = phenotype_datasets,
     target = "rules", 
     rules_selected = "best", 
     cut_off = c(0.5, 0, 0.5, 0), 
     measures_used = c("cosine", "jaccard", "kulczynski", "support"), 
     data_partition = data_partition,
     title = "Excluding Streptomycin"
 )


indicator_2017_comparison_plot <- gene_vs_phenotype_cumulative_percent_captured_multi_v3(
     datasets_geno = genotype2017_datasets,
     datasets_pheno = phenotype2017_datasets,
     target = "rules", 
     rules_selected = "best", 
     cut_off = c(0.5, 0, 0.5, 0), 
     measures_used = c("cosine", "jaccard", "kulczynski", "support"), 
     data_partition = data_partition,
     title = "Including Streptomycin"
 )




indicator_comparison_plot <- format_and_combine_plots(
    plot_list = list(indicator_comparison_plot, indicator_2017_comparison_plot),
    title_text = "",  # Remove title
    labels = c("E", "F"),
    ncol = 2,
    nrow = 1
)
print("DONE WITH BETWEEN INDICATORS PLOT")

#####################################
#RUN THE BETWEEN DATASETS COMPARISON
#####################################
print("STARTED BETWEEN SOURCES PLOT")

#only keep years in phenotype dfs that are present in both cecal and retail
cecal_phenotype_df <- cecal_phenotype_df[cecal_phenotype_df$Year %in% Retail_Meats_phenotype_df$Year, ]
Retail_Meats_phenotype_df <- Retail_Meats_phenotype_df[Retail_Meats_phenotype_df$Year %in% cecal_phenotype_df$Year, ]
#do this for the class level dfs too
cecal_Class_phenotype_df <- cecal_Class_phenotype_df[cecal_Class_phenotype_df$Year %in% Retail_Meats_Class_phenotype_df$Year, ]
Retail_Meats_Class_phenotype_df <- Retail_Meats_Class_phenotype_df[Retail_Meats_Class_phenotype_df$Year %in% cecal_Class_phenotype_df$Year, ]


#only keep years in genotype dfs that are present in both cecal and retail
cecal_Family_genotype_df <- cecal_Family_genotype_df[cecal_Family_genotype_df$Year %in% Retail_Meats_Family_genotype_df$Year, ]
Retail_Meats_Family_genotype_df <- Retail_Meats_Family_genotype_df[Retail_Meats_Family_genotype_df$Year %in% cecal_Family_genotype_df$Year, ]
#do this for the class level dfs too
cecal_Class_genotype_df <- cecal_Class_genotype_df[cecal_Class_genotype_df$Year %in% Retail_Meats_Class_genotype_df$Year, ]
Retail_Meats_Class_genotype_df <- Retail_Meats_Class_genotype_df[Retail_Meats_Class_genotype_df$Year %in% cecal_Class_genotype_df$Year, ]


phenotype_datasets <- list(
    "Cecal Phenotype" = cecal_Class_phenotype_df,
    "Retail Meats Phenotype" = Retail_Meats_Class_phenotype_df,
    "NAHLN Phenotype" = NAHLN_Class_phenotype_df
)

between_sources_phenotype_plot <- phenotype_cumulative_percent_captured_multi(
    datasets_pheno = phenotype_datasets,
    target = "rules", 
    rules_selected = "best", 
    cut_off = c(0.5, 0, 0.5, 0), 
    measures_used = c("cosine", "jaccard", "kulczynski", "support"), 
    data_partition = data_partition,
    class_level = TRUE
)





# Create the dataset lists
genotype_datasets <- list(
    "Cecal Genotype" = cecal_Class_genotype_df,
    "Retail Meats Genotype" = Retail_Meats_Class_genotype_df,
    "NAHLN Genotype" = NAHLN_Class_genotype_df
)





between_sources_genotype_plot <- genotype_cumulative_percent_captured_multi(
    datasets_geno = genotype_datasets,
    target = "rules", 
    rules_selected = "best", 
    cut_off = c(0.5, 0, 0.5, 0), 
    measures_used = c("cosine", "jaccard", "kulczynski", "support"), 
    data_partition = data_partition,
    class_level = TRUE
)


between_data_sources_plot <- format_and_combine_plots(
    plot_list = list(between_sources_phenotype_plot, between_sources_genotype_plot),
    title_text = "",  # Remove title
    labels = c("C", "D"),
    legend_size = 10
)
print("DONE WITH BETWEEN SOURCES PLOT")

# Save the final plot with increased dimensions
# ggsave(plot = between_data_sources_plot, filename = str_glue("combined_figures_outputs/percent_captured_analysis/between_data_sources/rule_set_comparison_between_data_sources.png"), 
#        width = 18, height = 10, dpi = 300)



#################################
####### MAKE FINAL PANEL ########
#################################


#stack three plots (within_indicators_plot, between_data_sources_plot, and indicator_comparison_plot) on top of each other to make one big panel

# Create the final combined panel of all plots
final_panel <- plot_grid(
    within_indicators_plot,
    between_data_sources_plot,
    indicator_comparison_plot,
    ncol = 1,
    nrow = 3,
    rel_heights = c(1, 1, 1),
    align = "v",
    axis = "lr"
)

# Add border to final panel with padding (no title)
final_panel <- ggdraw(final_panel) +
    theme(
        plot.background = element_rect(color = "black", fill = NA, linewidth = 3),
        plot.margin = margin(20, 20, 20, 20)
    )

# Save with larger dimensions to accommodate the border
ggsave(
    plot = final_panel,
    filename = str_glue("combined_outputs/combined_figures_outputs/combined_analysis_panel_Fig2_in_manuscript.png"),
    width = 18,
    height = 30,
    dpi = 300,
    bg = "white"  # Ensure white background
)




}



































#####################################################################################################################################################
#####################################################################################################################################################
######################################### RUN ANALYSES  #############################################################################################
#####################################################################################################################################################
#####################################################################################################################################################

#load in the data

#cecal and retail meats
cecal_phenotype_df <- read.csv("dataset_specific_outputs/cecal/NEW_cecal_wide_resStatus_phenotype.csv")
cecal_Class_phenotype_df <- read.csv("dataset_specific_outputs/cecal/NEW_cecal_wide_class_level_phenotype.csv")

cecal_Family_genotype_df <- read.csv("dataset_specific_outputs/cecal/cecal_wide_corrected_genotype.csv")
cecal_Class_genotype_df <- read.csv("dataset_specific_outputs/cecal/cecal_wide_class_level_genotype.csv")

Retail_Meats_phenotype_df <- read.csv("dataset_specific_outputs/Retail_Meats/NEW_Retail_Meats_wide_resStatus_phenotype.csv")
Retail_Meats_Class_phenotype_df <- read.csv("dataset_specific_outputs/Retail_Meats/NEW_Retail_Meats_wide_class_level_phenotype.csv")

Retail_Meats_Family_genotype_df <- read.csv("dataset_specific_outputs/Retail_Meats/Retail_Meats_wide_corrected_genotype.csv")
Retail_Meats_Class_genotype_df <- read.csv("dataset_specific_outputs/Retail_Meats/Retail_Meats_wide_class_level_genotype.csv")

#cecal and retail meats 2017 onward
cecal_2017_phenotype_df <- read.csv("dataset_specific_outputs/cecal_2017Onward/NEW_cecal_2017Onward_wide_resStatus_phenotype.csv")
cecal_2017_Class_phenotype_df <- read.csv("dataset_specific_outputs/cecal_2017Onward/NEW_cecal_2017Onward_wide_class_level_phenotype.csv")

cecal_2017_Family_genotype_df <- read.csv("dataset_specific_outputs/cecal_2017Onward/cecal_2017Onward_wide_corrected_genotype.csv")
cecal_2017_Class_genotype_df <- read.csv("dataset_specific_outputs/cecal_2017Onward/cecal_2017Onward_wide_class_level_genotype.csv")

Retail_2017_Meats_phenotype_df <- read.csv("dataset_specific_outputs/Retail_Meats_2017Onward/NEW_Retail_Meats_2017Onward_wide_resStatus_phenotype.csv")
Retail_2017_Meats_Class_phenotype_df <- read.csv("dataset_specific_outputs/Retail_Meats_2017Onward/NEW_Retail_Meats_2017Onward_wide_class_level_phenotype.csv")

Retail_2017_Meats_Family_genotype_df <- read.csv("dataset_specific_outputs/Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_corrected_genotype.csv")
Retail_2017_Meats_Class_genotype_df <- read.csv("dataset_specific_outputs/Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_class_level_genotype.csv")

#NAHLN
NAHLN_phenotype_df <- read.csv("dataset_specific_outputs/NAHLN/NAHLN_wide_resStatus_phenotype.csv")
NAHLN_Class_phenotype_df <- read.csv("dataset_specific_outputs/NAHLN/NAHLN_wide_class_level_phenotype.csv")

NAHLN_Family_genotype_df <- read.csv("dataset_specific_outputs/NAHLN/NAHLN_wide_corrected_genotype.csv")
NAHLN_Class_genotype_df <-read.csv("dataset_specific_outputs/NAHLN/NAHLN_wide_class_level_genotype.csv")


run_all_percent_captured_analyses(data_partition = "full_date_range",

                                            cecal_phenotype_df, 
                                            cecal_Class_phenotype_df,

                                            cecal_Family_genotype_df ,
                                            cecal_Class_genotype_df,

                                            Retail_Meats_phenotype_df,
                                            Retail_Meats_Class_phenotype_df,

                                            Retail_Meats_Family_genotype_df,
                                            Retail_Meats_Class_genotype_df,


                                            cecal_2017_phenotype_df, 
                                            cecal_2017_Class_phenotype_df,

                                            cecal_2017_Family_genotype_df ,
                                            cecal_2017_Class_genotype_df,

                                            Retail_2017_Meats_phenotype_df,
                                            Retail_2017_Meats_Class_phenotype_df,

                                            Retail_2017_Meats_Family_genotype_df,
                                            Retail_2017_Meats_Class_genotype_df,


                                            
                                            NAHLN_phenotype_df, 
                                            NAHLN_Class_phenotype_df,

                                            NAHLN_Family_genotype_df ,
                                            NAHLN_Class_genotype_df


                                            )



print("DONE RUNNING THE COMBINED PERCENT CAPTURED ANALYSIS")