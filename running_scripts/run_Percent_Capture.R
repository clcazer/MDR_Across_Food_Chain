
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
                                            Retail_Meats_Class_genotype_df
                                                ) {

    #####################################
#RUN THE WITHIN INDICATORS COMPARISON
#####################################
phenotype_datasets <- list(
    "Cecal Phenotype" = cecal_phenotype_df,
    "Retail Meats Phenotype" = Retail_Meats_phenotype_df
)

class_phenotype_datasets <- list(
    "Cecal Class Phenotype" = cecal_Class_phenotype_df,
    "Retail Meats Class Phenotype" = Retail_Meats_Class_phenotype_df
)

genotype_datasets <- list(
    "Cecal Genotype" = cecal_Family_genotype_df,
    "Retail Meats Genotype" = Retail_Meats_Family_genotype_df
)

class_genotype_datasets <- list(
    "Cecal Class Genotype" = cecal_Class_genotype_df,
    "Retail Meats Class Genotype" = Retail_Meats_Class_genotype_df
)

 nonclass_plot <- cumulative_proportion_captured_multi(
    datasets_geno = genotype_datasets,
    datasets_pheno = phenotype_datasets,
    class_level = FALSE,
    target = "rules",
    rules_selected = "best",
    cut_off = c(0.5, 0, 0.5, 0),
    measures_used = c("cosine", "jaccard", "kulczynski", "support"),
    data_partition = data_partition
)

 class_plot <- cumulative_proportion_captured_multi(
    datasets_geno = class_genotype_datasets,
    datasets_pheno = class_phenotype_datasets,
    class_level = TRUE,
    target = "rules",
    rules_selected = "best",
    cut_off = c(0.5, 0, 0.5, 0),
    measures_used = c("cosine", "jaccard", "kulczynski", "support"),
    data_partition = data_partition
)






format_and_combine_plots <- function(plot_list, title_text, labels, ncol = 2, nrow =1, title = TRUE) {
    # Format individual plots
    formatted_plots <- lapply(seq_along(plot_list), function(i) {
        p <- plot_list[[i]]
        p <- p + 
            theme(
                plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 10)),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
                axis.title.x = element_text(margin = margin(t = 10)),
                axis.title.y = element_text(margin = margin(r = 10)),
                #plot.background = element_rect(color = "black", fill = NA), # Add border to individual plots
                #panel.border = element_rect(color = "black", fill = NA)
            ) +
            scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1))
        
        # Remove y-axis title for right column plots
        if (i %% 2 == 0) {
            p <- p + theme(axis.title.y = element_blank())
        }
        
        return(p)
    })
    
    # Create panel of plots
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

    # Add border to combined plot
   # combined_plot <- ggdraw(combined_plot) + 
      #  draw_plot_label(label = "", size = 11, fontface = "bold") +
     #   theme(
      #  plot.background = element_rect(color = "black", fill = NA, linewidth = 1),
       # plot.margin = margin(20, 20, 20, 20)
       # )

    if (title == TRUE) {
        # Add title
        title <- ggdraw() + 
            draw_label(title_text,
                      fontface = "bold",
                      size = 25,
                      colour = "black")
        
        # Combine title and plots
        final_plot <- plot_grid(
            title, combined_plot,
            ncol = 1,
            rel_heights = c(0.1, 1)
        )
        
        # Add border to final plot
       # final_plot <- ggdraw(final_plot) +
         #   theme(
       # plot.background = element_rect(color = "black", fill = NA, linewidth = 1),
        #plot.margin = margin(20, 20, 20, 20)  # Add padding around the plot
        #)
        
        return(final_plot)
    } else {
        return(combined_plot)
    }
}

within_indicators_plot <- format_and_combine_plots(
    plot_list = list(nonclass_plot, class_plot),
    title_text = "Rule Set Comparison Within Indicators",
    labels = c("A", "B")
)

# Save the final plot with increased dimensions
ggsave(plot = within_indicators_plot, filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/within_indicators/rule_set_comparison_within_indicators.png"), 
       width = 18, height = 10, dpi = 300)


#####################################
#RUN THE BETWEEN INDICATORS COMPARISON
#####################################
#since we're comparing between indicators, data has to be at the class level

#Include only rows where the ID matches in both dataframes for each data source
IDMathced_Retail_Meats_Class_phenotype_df <- Retail_Meats_Class_phenotype_df[Retail_Meats_Class_phenotype_df$ID %in% Retail_Meats_Class_genotype_df$ID, ]
IDMathced_Retail_Meats_Class_genotype_df <- Retail_Meats_Class_genotype_df[Retail_Meats_Class_genotype_df$ID %in% IDMathced_Retail_Meats_Class_phenotype_df$ID, ]

IDMathced_cecal_Class_phenotype_df <- cecal_Class_phenotype_df[cecal_Class_phenotype_df$ID %in% cecal_Class_genotype_df$ID, ]
IDMathced_cecal_Class_genotype_df <- cecal_Class_genotype_df[cecal_Class_genotype_df$ID %in% IDMathced_cecal_Class_phenotype_df$ID, ]

indicator_comparison_plot <- gene_vs_phenotype_cumulative_percent_captured_multi(
     genotype_df1 = IDMathced_Retail_Meats_Class_genotype_df, 
     phenotype_df1 = IDMathced_Retail_Meats_Class_phenotype_df, 
     genotype_df2 = IDMathced_cecal_Class_genotype_df, 
     phenotype_df2 = IDMathced_cecal_Class_phenotype_df, 
     target = "rules", 
     rules_selected = "best", 
     cut_off = c(0.5, 0, 0.5, 0), 
     measures_used = c("cosine", "jaccard", "kulczynski", "support"), 
     data_source1 = "Retail_Meats",
     data_source2 = "Cecal",
     data_partition = data_partition
 )


indicator_comparison_plot <- format_and_combine_plots(
    plot_list = list(indicator_comparison_plot),
    title_text = "Rule Set Comparison Between Indicators",
    labels = c("E"),
    ncol = 1,
    nrow = 1,
    #title = FALSE
)

#####################################
#RUN THE BETWEEN DATASETS COMPARISON
#####################################
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



data_source_nonclass_plot <- cecal_vs_retail_cumulative_percent_captured_multi(
    cecal_genotype_df = cecal_Family_genotype_df, 
    cecal_phenotype_df = cecal_phenotype_df, 
    retail_genotype_df = Retail_Meats_Family_genotype_df, 
    retail_phenotype_df = Retail_Meats_phenotype_df, 
    target = "rules", 
    rules_selected = "best", 
    cut_off = c(0.5, 0, 0.5, 0), 
    measures_used = c("cosine", "jaccard", "kulczynski", "support"), 
    data_source1 = "Cecal",
    data_source2 = "Retail Meats",
    data_partition = data_partition,
    class_level = FALSE
)


data_source_class_plot <- cecal_vs_retail_cumulative_percent_captured_multi(
    cecal_genotype_df = cecal_Class_genotype_df, 
    cecal_phenotype_df = cecal_Class_phenotype_df, 
    retail_genotype_df = Retail_Meats_Class_genotype_df, 
    retail_phenotype_df = Retail_Meats_Class_phenotype_df, 
    target = "rules", 
    rules_selected = "best", 
    cut_off = c(0.5, 0, 0.5, 0), 
    measures_used = c("cosine", "jaccard", "kulczynski", "support"), 
    data_source1 = "Cecal",
    data_source2 = "Retail Meats",
    data_partition = data_partition,
    class_level = TRUE
)


between_data_sources_plot <- format_and_combine_plots(
    plot_list = list(data_source_nonclass_plot, data_source_class_plot),
    title_text = "Rule Set Comparison Between Sources",
    labels = c("C", "D")
)

# Save the final plot with increased dimensions
ggsave(plot = between_data_sources_plot, filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/between_data_sources/rule_set_comparison_between_data_sources.png"), 
       width = 18, height = 10, dpi = 300)




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

# Add border to final panel with padding
final_panel <- ggdraw(final_panel) +
    theme(
        plot.background = element_rect(color = "black", fill = NA, linewidth = 3),
        plot.margin = margin(20, 20, 20, 20)  # Add padding around the plot
    )

# Save with larger dimensions to accommodate the border
ggsave(
    plot = final_panel,
    filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/combined_analysis_panel.png"),
    width = 18,
    height = 30,
    dpi = 300,
    bg = "white"  # Ensure white background
)




}



#####################################################################################################################################################
#####################################################################################################################################################
######################################### RUN ANALYSES FOR THE FULL DATE RANGE ######################################################################
#####################################################################################################################################################
#####################################################################################################################################################

#load in the data
cecal_phenotype_df <- read.csv("cecal/cecal_wide_resStatus_phenotype.csv")
cecal_Class_phenotype_df <- read.csv("cecal/cecal_wide_class_level_phenotype.csv")

cecal_Family_genotype_df <- read.csv("cecal/cecal_wide_family_level_genotype.csv")
cecal_Class_genotype_df <- read.csv("cecal/cecal_wide_class_level_genotype.csv")

Retail_Meats_phenotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_resStatus_phenotype.csv")
Retail_Meats_Class_phenotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_phenotype.csv")

Retail_Meats_Family_genotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_family_level_genotype.csv")
Retail_Meats_Class_genotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_genotype.csv")


run_all_percent_captured_analyses(data_partition = "full_date_range",

                                            cecal_phenotype_df, 
                                            cecal_Class_phenotype_df,

                                            cecal_Family_genotype_df ,
                                            cecal_Class_genotype_df,

                                            Retail_Meats_phenotype_df,
                                            Retail_Meats_Class_phenotype_df,

                                            Retail_Meats_Family_genotype_df,
                                            Retail_Meats_Class_genotype_df
                                            )






#####################################################################################################################################################
#####################################################################################################################################################
######################################### RUN ANALYSES FOR THE FULL 2017 ONWARD #####################################################################
#####################################################################################################################################################
#####################################################################################################################################################



#load in the data
cecal_phenotype_df2017 <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_resStatus_phenotype.csv")
cecal_Class_phenotype_df2017 <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_class_level_phenotype.csv")

cecal_Family_genotype_df2017 <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_family_level_genotype.csv")
cecal_Class_genotype_df2017 <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_class_level_genotype.csv")

Retail_Meats_phenotype_df2017 <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_resStatus_phenotype.csv")
Retail_Meats_Class_phenotype_df2017 <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_class_level_phenotype.csv")

Retail_Meats_Family_genotype_df2017 <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_family_level_genotype.csv")
Retail_Meats_Class_genotype_df2017 <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_class_level_genotype.csv")


run_all_percent_captured_analyses(data_partition = "2017_onward",

                                            cecal_phenotype_df = cecal_phenotype_df2017, 
                                            cecal_Class_phenotype_df = cecal_Class_phenotype_df2017,

                                            cecal_Family_genotype_df = cecal_Family_genotype_df2017,
                                            cecal_Class_genotype_df = cecal_Class_genotype_df2017,

                                            Retail_Meats_phenotype_df = Retail_Meats_phenotype_df2017,
                                            Retail_Meats_Class_phenotype_df = Retail_Meats_Class_phenotype_df2017,

                                            Retail_Meats_Family_genotype_df = Retail_Meats_Family_genotype_df2017,
                                            Retail_Meats_Class_genotype_df = Retail_Meats_Class_genotype_df2017
                                            )











#####################################################################################################################################################
#####################################################################################################################################################
###################################### RUN ANALYSES FOR THE SOURCE ENCODINGS ADDED ##################################################################
#####################################################################################################################################################
#####################################################################################################################################################


