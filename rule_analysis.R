# Load necessary libraries
library(arules)
library(stringr)
library(rlist)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(readxl)
library(cowplot)
library(arulesViz)

# Source additional R scripts
source("data_wrangling.R")
source("rule_mining_and_selection.R")

# Function to create scatter plot comparing number of rules before and after implementing measure cut-offs
plot_all_vs_best_rules <- function(df, resistance_indicator, target, cut_off, measures_used, data_source) {
    # Initialize lists to store data
    all_rules_length <- list()
    best_rules_length <- list()
    year_list <- list()
    sample_size <- list()

    # Iterate through each year in the dataset
    for (year in c(min(df[, "Year"]):max(df[, "Year"]))) {
        # Get all rules for the current year
        all_rules <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator,
                                           year = year, target = target)
        
        all_rules_length <- list.append(all_rules_length, length(all_rules))

        # Apply cut-offs to get best rules
        best_rules <- implement_cuts(df = df, resistance_indicator = resistance_indicator, 
                                     target = target, cut_off = cut_off, year = year, measures_used = measures_used)
        
        best_rules_length <- list.append(best_rules_length, length(best_rules))
        year_list <- list.append(year_list, year)

        # Calculate sample size for the current year
        sample_size_df <- df[df$Year == year,]
        sample_size <- list.append(sample_size, length(unique(sample_size_df$ID)))
    }

    # Create dataframe with collected data
    best_v_all_df <- data.frame(Sample_Size = unlist(sample_size), Year = as.factor(unlist(year_list)), 
                                All_Rules = unlist(all_rules_length), Best_Rules = unlist(best_rules_length))

    print(best_v_all_df)

    # Create and save the plot
    ggplot(best_v_all_df, aes(x = All_Rules, y = Best_Rules, col = Year)) + 
        geom_point(size = 7) +
        theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black')) +
        labs(y = str_glue("Best Rules"), x = str_glue("All Rules"), title = str_glue("Rule Pruning Results"))
    
    ggsave(str_glue("{data_source}/figures/pruning_results/{resistance_indicator}/{data_source}_{resistance_indicator}_{target}_Rule_Pruning_Results.png"))
}








################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
#####################################################   WITHIN INDICATORS      ################################################################################################
################################################################################################################################################################################    
################################################################################################################################################################################










# Function to compare rule overlap between two consecutive years
rule_overlap <- function(df, resistance_indicator, class_level = FALSE, target, rules_selected = "best", cut_off, measures_used, data_source) {
    overlap_list <- list()
    year_list <- list()

    # Iterate through years, comparing each year with the next
    for (year in c(min(df[, "Year"]):(max(df[, "Year"])-1))) {
        # Get rules for current and next year based on selection criteria
        if (rules_selected == "best") {
            rules1 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = year)
            rules2 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = (year + 1))
        } else if (rules_selected == "all") {
            rules1 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = year, target = target)
            rules2 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = (year + 1), target = target)
        }

        # Fix item encodings to make them comparable across years
        new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
        if (target == "rules") {
            rules1_fixed <- new("rules",
                lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
            )
            rules2_fixed <- new("rules",
                lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
            )
        } else {
            rules1_fixed <- rules1
            rules2_fixed <- rules2
        }

        # Calculate overlap and store results
        overlap <- length((intersect(rules1_fixed, rules2_fixed)))/length((union(rules1_fixed, rules2_fixed)))
        overlap_list <- list.append(overlap_list, overlap)
        year_list <- list.append(year_list, year)
    }

    # Create dataframe with overlap data
    overlap_df <- data.frame(rules_shared = unlist(overlap_list), year = unlist(year_list))

    # Create and save the plot
    plot1 <- ggplot(overlap_df, aes(x = year, y = rules_shared)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black')) +
        labs(y = str_glue("Proportion Rules Shared"), x = "Year", title = str_glue("{resistance_indicator} Rules Shared")) +
        scale_x_continuous(breaks = round(seq(min(overlap_df$year), max(overlap_df$year), by = 2),1)) +
        ylim(0,1)
    if (class_level == FALSE) {
    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/{resistance_indicator}/proportion_shared/{data_source}_{resistance_indicator}_{target}_prop_shared_{rules_selected}.png"))
    } else {
    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/{resistance_indicator}/proportion_shared/{data_source}_{resistance_indicator}_{target}_prop_shared_CLASS_{rules_selected}.png"))
    }
}







# Function to compare rule overlap between two consecutive years
proportion_captured <- function(df, resistance_indicator, class_level = FALSE, target, rules_selected = "best", cut_off, measures_used, data_source) {
    overlap_list <- list()
    year_list <- list()

    # Iterate through years, comparing each year with the next
    for (year in c(min(df[, "Year"]):(max(df[, "Year"])-1))) {
        # Get rules for current and next year based on selection criteria
        if (rules_selected == "best") {
            rules1 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = year)
            rules2 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = (year + 1))
        } else if (rules_selected == "all") {
            rules1 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = year, target = target)
            rules2 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = (year + 1), target = target)
        }

        # Fix item encodings to make them comparable across years
        new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
        if (target == "rules") {
            rules1_fixed <- new("rules",
                lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
            )
            rules2_fixed <- new("rules",
                lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
            )
        } else {
            rules1_fixed <- rules1
            rules2_fixed <- rules2
        }

        # Calculate overlap and store results
        overlap <- length((intersect(rules1_fixed, rules2_fixed)))/length(((rules1_fixed)))
        overlap_list <- list.append(overlap_list, overlap)
        year_list <- list.append(year_list, year)
    }

    # Create dataframe with overlap data
    overlap_df <- data.frame(rules_shared = unlist(overlap_list), year = unlist(year_list))

    # Create and save the plot
    plot1 <- ggplot(overlap_df, aes(x = year, y = rules_shared)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black')) +
        labs(y = str_glue("Proportion Rules Captured in Subsequent Year"), x = "Year", title = str_glue("{resistance_indicator} Rules Shared")) +
        scale_x_continuous(breaks = round(seq(min(overlap_df$year), max(overlap_df$year), by = 2),1)) +
        ylim(0,1)
    if (class_level == FALSE) {
    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/{resistance_indicator}/proportion_captured/{data_source}_{resistance_indicator}_{target}_prop_captured_{rules_selected}.png"))
    } else {
    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/{resistance_indicator}/proportion_captured/{data_source}_{resistance_indicator}_{target}_prop_captured_CLASS_{rules_selected}.png"))
    }
}













# Function to calculate and plot cumulative average of proportion rules shared
cumulative_rule_stability <- function(df, resistance_indicator,class_level = FALSE, target, rules_selected = "best", cut_off, measures_used, data_source) {
    crs_list <- list()
    year_list <- list()

    for (year in c(min(df[, "Year"]):(max(df[, "Year"])-1))) {
        i <- match(year, sort(unique(df$Year)))

        if (rules_selected == "best") {
            rules1 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = year)
            rules2 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = (year + 1))
        } else if (rules_selected == "all") {
            rules1 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = year, target = target)
            rules2 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = (year + 1), target = target)
        }

        # Fix item encodings to make them comparable across years
        new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
        if (target == "rules") {
            rules1_fixed <- new("rules",
                lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
            )
            rules2_fixed <- new("rules",
                lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
            )
        } else {
            rules1_fixed <- rules1
            rules2_fixed <- rules2
        }

        val1 <- length(union(rules1_fixed, rules2_fixed))
        val2 <- length(rules1_fixed) + length(rules2_fixed) - length(intersect(rules1_fixed, rules2_fixed))
        print(val1 == val2)
        print(all.equal(itemLabels(rules1_fixed), itemLabels(rules2_fixed))) 
        print("======")  
    
        if (i == 1) {
            crs <- length((intersect(rules1_fixed, rules2_fixed)))/length((union(rules1_fixed, rules2_fixed)))
        } else {
            crs <- (1/i)*((i-1)*crs + length((intersect(rules1_fixed, rules2_fixed)))/length((union(rules1_fixed, rules2_fixed))))
        }
        print(crs)
        crs_list <- list.append(crs_list, crs)
        year_list <- list.append(year_list, year)
    }

    crs_df <- data.frame(crs = unlist(crs_list), year = unlist(year_list))

    plot1 <- ggplot(crs_df, aes(x = year, y = crs)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black')) +
        labs(y = str_glue("Cumulative Rule Stability"), x = "Year", title = str_glue("{resistance_indicator} crs")) +
        scale_x_continuous(breaks = round(seq(min(crs_df$year), max(crs_df$year), by = 2),1)) +
        ylim(0,1)
    if (class_level == FALSE) {
    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/{resistance_indicator}/crs/{data_source}_{resistance_indicator}_{target}_crs_{rules_selected}.png"))
    } else {
    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/{resistance_indicator}/crs/{data_source}_{resistance_indicator}_{target}_crs_CLASS_{rules_selected}.png"))
    }
}





# Function to calculate and plot cumulative average of proportion rules shared
cumulative_proportion_captured <- function(df, resistance_indicator, class_level = FALSE, target, rules_selected = "best", cut_off, measures_used, data_source) {
        crs_list <- list()
        year_list <- list()

        for (year in c(min(df[, "Year"]):(max(df[, "Year"])-1))) {
            i <- match(year, sort(unique(df$Year)))

            if (rules_selected == "best") {
                rules1 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = year)
                rules2 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = (year + 1))
            } else if (rules_selected == "all") {
                rules1 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = year, target = target)
                rules2 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = (year + 1), target = target)
            }

            # Fix item encodings to make them comparable across years
            new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
            if (target == "rules") {
                rules1_fixed <- new("rules",
                    lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                    rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
                )
                rules2_fixed <- new("rules",
                    lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                    rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
                )
            } else {
                rules1_fixed <- rules1
                rules2_fixed <- rules2
            }

        val1 <- length(union(rules1_fixed, rules2_fixed))
        val2 <- length(rules1_fixed) + length(rules2_fixed) - length(intersect(rules1_fixed, rules2_fixed))
        print(val1 == val2)
        print(all.equal(itemLabels(rules1_fixed), itemLabels(rules2_fixed))) 
        print("======")  

            if (i == 1) {
                crs <- length((intersect(rules1_fixed, rules2_fixed)))/length((rules1_fixed))
            } else {
                crs <- (1/i)*((i-1)*crs + length((intersect(rules1_fixed, rules2_fixed)))/length(rules1_fixed))
            }
        print(crs)
            crs_list <- list.append(crs_list, crs)
            year_list <- list.append(year_list, year)
        }

    crs_df <- data.frame(crs = unlist(crs_list), year = unlist(year_list))

    plot1 <- ggplot(crs_df, aes(x = year, y = crs)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black')) +
        labs(y = str_glue("Cumulative Proportion Captured"), x = "Year", title = str_glue("{resistance_indicator} crs")) +
        scale_x_continuous(breaks = round(seq(min(crs_df$year), max(crs_df$year), by = 2),1)) +
        ylim(0,1)
    if (class_level == FALSE) {
    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/{resistance_indicator}/c_proportion_captured/{data_source}_{resistance_indicator}_{target}_c_prop_captured_{rules_selected}.png"))
    } else {
    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/{resistance_indicator}/c_proportion_captured/{data_source}_{resistance_indicator}_{target}_c_prop_captured_CLASS_{rules_selected}.png"))
    }
}

# Function to calculate and plot cumulative average of proportion rules shared for multiple datasets
cumulative_proportion_captured_multi <- function(datasets, resistance_indicator, class_level = FALSE, target, rules_selected = "best", cut_off, measures_used, data_source) {
    all_crs_df <- data.frame()

    for (dataset_name in names(datasets)) {
        df <- datasets[[dataset_name]]
        crs_list <- list()
        year_list <- list()

        for (year in c(min(df[, "Year"]):(max(df[, "Year"])-1))) {
            i <- match(year, sort(unique(df$Year)))

            if (rules_selected == "best") {
                rules1 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = year)
                rules2 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = (year + 1))
            } else if (rules_selected == "all") {
                rules1 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = year, target = target)
                rules2 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = (year + 1), target = target)
            }

            # Fix item encodings to make them comparable across years
            new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
            if (target == "rules") {
                rules1_fixed <- new("rules",
                    lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                    rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
                )
                rules2_fixed <- new("rules",
                    lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                    rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
                )
            } else {
                rules1_fixed <- rules1
                rules2_fixed <- rules2
            }

            if (i == 1) {
                crs <- length((intersect(rules1_fixed, rules2_fixed)))/length((rules1_fixed))
            } else {
                crs <- (1/i)*((i-1)*crs + length((intersect(rules1_fixed, rules2_fixed)))/length(rules1_fixed))
            }
            crs_list <- list.append(crs_list, crs)
            year_list <- list.append(year_list, year)
        }

        crs_df <- data.frame(crs = unlist(crs_list), year = unlist(year_list), dataset = dataset_name)
        all_crs_df <- rbind(all_crs_df, crs_df)
    }
if (resistance_indicator == "phenotype") {
    plot1 <- ggplot(all_crs_df, aes(x = year, y = crs, color = dataset)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black'),
              legend.position = "bottom") +
        labs(y = str_glue("Cumulative Proportion Captured"), x = "Year", 
             title = str_glue("Phenotype Cumulative Proportion Captured"),
             color = "Dataset") +
        scale_x_continuous(breaks = round(seq(min(all_crs_df$year), max(all_crs_df$year), by = 2),1)) +
        ylim(0,1)
} else {
plot1 <- ggplot(all_crs_df, aes(x = year, y = crs, color = dataset)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black'),
              legend.position = "bottom") +
        labs(y = str_glue("Cumulative Proportion Captured"), x = "Year", 
             title = str_glue("Genotype Cumulative Proportion Captured"),
             color = "Dataset") +
        scale_x_continuous(breaks = round(seq(min(all_crs_df$year), max(all_crs_df$year), by = 2),1)) +
        ylim(0,1)
}
    

    filename_suffix <- if(class_level) "CLASS" else ""
    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/{resistance_indicator}/c_proportion_captured/multi_{resistance_indicator}_{target}_c_prop_captured_{rules_selected}_{filename_suffix}.png"))
    
    return(plot1)
}

cecal_phenotype_df <- read.csv("cecal/cecal_wide_resStatus_phenotype.csv")
cecal_Class_phenotype_df <- read.csv("cecal/cecal_wide_class_level_phenotype.csv")
cecal_Family_genotype_df <- read.csv("cecal/cecal_wide_family_level_genotype.csv")
cecal_Class_genotype_df <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_class_level_genotype.csv")

Retail_Meats_phenotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_resStatus_phenotype.csv")
Retail_Meats_Class_phenotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_phenotype.csv")
Retail_Meats_Family_genotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_family_level_genotype.csv")
Retail_Meats_Class_genotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_genotype.csv")

#filter all dataframes to only include years 2013 and onward
cecal_phenotype_df <- cecal_phenotype_df[cecal_phenotype_df$Year >= 2013,]
cecal_Class_phenotype_df <- cecal_Class_phenotype_df[cecal_Class_phenotype_df$Year >= 2013,]
cecal_Family_genotype_df <- cecal_Family_genotype_df[cecal_Family_genotype_df$Year >= 2013,]
cecal_Class_genotype_df <- cecal_Class_genotype_df[cecal_Class_genotype_df$Year >= 2013,]

Retail_Meats_phenotype_df <- Retail_Meats_phenotype_df[Retail_Meats_phenotype_df$Year >= 2013,]
Retail_Meats_Class_phenotype_df <- Retail_Meats_Class_phenotype_df[Retail_Meats_Class_phenotype_df$Year >= 2013,]
Retail_Meats_Family_genotype_df <- Retail_Meats_Family_genotype_df[Retail_Meats_Family_genotype_df$Year >= 2013,]
Retail_Meats_Class_genotype_df <- Retail_Meats_Class_genotype_df[Retail_Meats_Class_genotype_df$Year >= 2013,]


datasets <- list(
    #"Cecal Class Genotype" = cecal_Class_genotype_df,
    #"Retail Meats Class Genotype" = Retail_Meats_Class_genotype_df
    "Cecal Class Phenotype" = cecal_phenotype_df,
    "Retail Meats Class Phenotype" = Retail_Meats_phenotype_df
)

plot <- cumulative_proportion_captured_multi(
    datasets = datasets,
    resistance_indicator = "phenotype",
    class_level = FALSE,
    target = "rules",
    rules_selected = "best",
    cut_off = c(0.5, 0, 0.5, 0),
    measures_used = c("cosine", "jaccard", "kulczynski", "support"),
    data_source = "all_combined"
)




# Function to compare rule overlap between genotype and phenotype
rule_overlap_gene_vs_phenotype <- function(genotype_df, phenotype_df, target, rules_selected = "best", cut_off, measures_used, data_source) {
    overlap_list <- list()
    year_list <- list()

    for (year in c(min(genotype_df[, "Year"]):(max(genotype_df[, "Year"])))) {
        if (rules_selected == "best") {
            rules1 <- implement_cuts(df = phenotype_df, resistance_indicator = "phenotype", target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = TRUE)
            rules2 <- implement_cuts(df = genotype_df, resistance_indicator = "genotype", target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = TRUE)
        } else if (rules_selected == "all") {
            rules1 <- get_rules_or_itemsets(df = phenotype_df, resistance_indicator = "phenotype", year = year, target = target, agg = TRUE)
            rules2 <- get_rules_or_itemsets(df = genotype_df, resistance_indicator = "genotype", year = year, target = target, agg = TRUE)
        }

        # Fix item encodings to make them comparable across years
        new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
        if (target == "rules") {
            rules1_fixed <- new("rules",
                lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
            )
            rules2_fixed <- new("rules",
                lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
            )
        } else {
            rules1_fixed <- rules1
            rules2_fixed <- rules2
        }

        overlap <- length((intersect(rules1_fixed, rules2_fixed)))/length((union(rules1_fixed, rules2_fixed)))
        overlap_list <- list.append(overlap_list, overlap)
        year_list <- list.append(year_list, year)
    }

    overlap_df <- data.frame(rules_shared = unlist(overlap_list), year = unlist(year_list))
    print(overlap_df)

    plot1 <- ggplot(overlap_df, aes(x = year, y = rules_shared)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black')) +
        labs(y = str_glue("Proportion Rules Shared"), x = "Year", title = str_glue("Genotype Vs. Phenotype Class-Level Rules Shared")) +
        scale_x_continuous(breaks = round(seq(min(overlap_df$year), max(overlap_df$year), by = 2),1)) +
        ylim(0,1)

    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/between_indicators/proportion_shared/{data_source}_{target}_GvP_prop_shared_{rules_selected}.png"))
}


















################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
#####################################################   BETWEEN INDICATORS      ################################################################################################
################################################################################################################################################################################    
################################################################################################################################################################################








# Function to calculate and plot cumulative rule stability between genotype and phenotype
cumulative_rule_stability_gene_vs_phenotype <- function(genotype_df, phenotype_df, target, rules_selected = "best", cut_off, measures_used, data_source) {
    crs_list <- list()
    year_list <- list()

    for (year in c(min(genotype_df[, "Year"]):(max(genotype_df[, "Year"])))) {
        i <- match(year, sort(unique(genotype_df$Year)))

        if (rules_selected == "best") {
            rules1 <- implement_cuts(df = phenotype_df, resistance_indicator = "phenotype", target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = TRUE)
            rules2 <- implement_cuts(df = genotype_df, resistance_indicator = "genotype", target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = TRUE)
        } else if (rules_selected == "all") {
            rules1 <- get_rules_or_itemsets(df = phenotype_df, resistance_indicator = "phenotype", year = year, target = target, agg = TRUE)
            rules2 <- get_rules_or_itemsets(df = genotype_df, resistance_indicator = "genotype", year = year, target = target, agg = TRUE)
        }

        # Fix item encodings to make them comparable across years
        new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
        if (target == "rules") {
            rules1_fixed <- new("rules",
                lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
            )
            rules2_fixed <- new("rules",
                lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
            )
        } else {
            rules1_fixed <- arules::recode(rules1, itemLabels = new_itemLabels)
            rules2_fixed <- arules::recode(rules2, itemLabels = new_itemLabels)
        }

        val1 <- length(union(rules1_fixed, rules2_fixed))
        val2 <- length(rules1_fixed) + length(rules2_fixed) - length(intersect(rules1_fixed, rules2_fixed))
        print(val1 == val2)
        print(all.equal(itemLabels(rules1_fixed), itemLabels(rules2_fixed))) 
        print("======") 
    
        if (i == 1) {
            crs <- length((intersect(rules1_fixed, rules2_fixed)))/length((union(rules1_fixed, rules2_fixed)))
        } else {
            crs <- (1/i)*((i-1)*crs + length((intersect(rules1_fixed, rules2_fixed)))/length((union(rules1_fixed, rules2_fixed))))
        }
        print(crs)
        crs_list <- list.append(crs_list, crs)
        year_list <- list.append(year_list, year)
    }

    crs_df <- data.frame(crs = unlist(crs_list), year = unlist(year_list))

    plot1 <- ggplot(crs_df, aes(x = year, y = crs)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black')) +
        labs(y = str_glue("Cumulative Rule Stability"), x = "Year", title = str_glue("Genotype Versus Phenotype crs")) +
        scale_x_continuous(breaks = round(seq(min(crs_df$year), max(crs_df$year), by = 2),1)) +
        ylim(0,1)

    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/between_indicators/crs/{data_source}_{target}_crs_{rules_selected}.png"))
}
























# Function to calculate and plot cumulative percent captured between genotype and phenotype
gene_vs_phenotype_cumulative_percent_captured <- function(genotype_df, phenotype_df, target, rules_selected = "best", cut_off, measures_used, data_source) {
    crs_list_geno <- list()
    crs_list_pheno <- list()
    year_list <- list()

    for (year in c(min(genotype_df[, "Year"]):(max(genotype_df[, "Year"])))) {
        i <- match(year, sort(unique(genotype_df$Year)))
        if (rules_selected == "best") {
            rules_pheno <- implement_cuts(df = phenotype_df, resistance_indicator = "phenotype", target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
            rules_geno <- implement_cuts(df = genotype_df, resistance_indicator = "genotype", target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
        } else if (rules_selected == "all") {
            rules_pheno <- get_rules_or_itemsets(df = phenotype_df, resistance_indicator = "phenotype", year = year, target = target, agg = FALSE)
            rules_geno <- get_rules_or_itemsets(df = genotype_df, resistance_indicator = "genotype", year = year, target = target, agg = FALSE)
        }

        # Fix item encodings to make them comparable across years
        new_itemLabels <- union(itemLabels(rules_pheno), itemLabels(rules_geno))
        if (target == "rules") {
            rules_pheno_fixed <- new("rules",
                lhs = arules::recode(lhs(rules_pheno), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules_pheno), itemLabels = new_itemLabels)
            )
            rules_geno_fixed <- new("rules",
                lhs = arules::recode(lhs(rules_geno), itemLabels = new_itemLabels), 
                rhs = arules::recode(rhs(rules_geno), itemLabels = new_itemLabels)
            )
        } else {
            rules_pheno_fixed <- arules::recode(rules_pheno, itemLabels = new_itemLabels)
            rules_geno_fixed <- arules::recode(rules_geno, itemLabels = new_itemLabels)
        }

        # Calculate percentages for both directions
        geno_in_pheno <- length(intersect(rules_geno_fixed, rules_pheno_fixed)) / length(rules_geno_fixed)
        pheno_in_geno <- length(intersect(rules_pheno_fixed, rules_geno_fixed)) / length(rules_pheno_fixed)

        if (i == 1) {
            crs_geno <- geno_in_pheno
            crs_pheno <- pheno_in_geno
        } else {
            crs_geno <- (1/i) * ((i-1) * crs_geno + geno_in_pheno)
            crs_pheno <- (1/i) * ((i-1) * crs_pheno + pheno_in_geno)
        }

        crs_list_geno <- list.append(crs_list_geno, crs_geno)
        crs_list_pheno <- list.append(crs_list_pheno, crs_pheno)
        year_list <- list.append(year_list, year)
    }

    crs_df <- data.frame(
        year = unlist(year_list),
        geno_in_pheno = unlist(crs_list_geno),
        pheno_in_geno = unlist(crs_list_pheno)
    )

    plot1 <- ggplot(crs_df) +
        geom_line(aes(x = year, y = geno_in_pheno, color = "Genotype in Phenotype"), linewidth = 2) +
        geom_point(aes(x = year, y = geno_in_pheno, color = "Genotype in Phenotype"), size = 3) +
        geom_line(aes(x = year, y = pheno_in_geno, color = "Phenotype in Genotype"), linewidth = 2) +
        geom_point(aes(x = year, y = pheno_in_geno, color = "Phenotype in Genotype"), size = 3) +
        theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black'),
              legend.position = "bottom") +
        labs(y = "Cumulative Proportion Captured", x = "Year", 
             title = "Cumulative Proportion of Rules Captured",
             color = "Comparison") +
        scale_x_continuous(breaks = round(seq(min(crs_df$year), max(crs_df$year), by = 2),1)) +
        scale_color_manual(values = c("Genotype in Phenotype" = "blue", "Phenotype in Genotype" = "red")) +
        ylim(0,1)

    ggsave(plot = plot1, filename = str_glue("{data_source}/figures/shared_rules_analysis/between_indicators/crs/{data_source}_{target}_captured_bidirectional_{rules_selected}.png"))
    
    return(plot1)
}


 genotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_genotype.csv")

#genotype_df <- subset(genotype_df, Year != 2020)

 phenotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_phenotype.csv")

# Include only rows where the ID matches in both dataframes
 phenotype_df <- phenotype_df[phenotype_df$ID %in% genotype_df$ID, ]
 genotype_df <- genotype_df[genotype_df$ID %in% phenotype_df$ID, ]


#gene_vs_phenotype_cumulative_percent_captured(genotype_df = genotype_df, phenotype_df = phenotype_df, target = "rules", rules_selected = "best", cut_off = c(0.5, 0, 0.5, 0), measures_used = c("cosine", "jaccard", "kulczynski", "support"), data_source = "Retail_Meats")








# Function to calculate and plot cumulative percent captured between genotype and phenotype for two data sources
gene_vs_phenotype_cumulative_percent_captured_multi <- function(genotype_df1, phenotype_df1, genotype_df2, phenotype_df2, target, rules_selected = "best", cut_off, measures_used, data_source1, data_source2) {
    process_dataset <- function(genotype_df, phenotype_df, data_source) {
        crs_list_geno <- list()
        crs_list_pheno <- list()
        year_list <- list()

        for (year in c(min(genotype_df[, "Year"]):(max(genotype_df[, "Year"])))) {
            i <- match(year, sort(unique(genotype_df$Year)))
            if (rules_selected == "best") {
                rules_pheno <- implement_cuts(df = phenotype_df, resistance_indicator = "phenotype", target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
                rules_geno <- implement_cuts(df = genotype_df, resistance_indicator = "genotype", target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
            } else if (rules_selected == "all") {
                rules_pheno <- get_rules_or_itemsets(df = phenotype_df, resistance_indicator = "phenotype", year = year, target = target, agg = FALSE)
                rules_geno <- get_rules_or_itemsets(df = genotype_df, resistance_indicator = "genotype", year = year, target = target, agg = FALSE)
            }

            new_itemLabels <- union(itemLabels(rules_pheno), itemLabels(rules_geno))
            if (target == "rules") {
                rules_pheno_fixed <- new("rules",
                    lhs = arules::recode(lhs(rules_pheno), itemLabels = new_itemLabels), 
                    rhs = arules::recode(rhs(rules_pheno), itemLabels = new_itemLabels)
                )
                rules_geno_fixed <- new("rules",
                    lhs = arules::recode(lhs(rules_geno), itemLabels = new_itemLabels), 
                    rhs = arules::recode(rhs(rules_geno), itemLabels = new_itemLabels)
                )
            } else {
                rules_pheno_fixed <- arules::recode(rules_pheno, itemLabels = new_itemLabels)
                rules_geno_fixed <- arules::recode(rules_geno, itemLabels = new_itemLabels)
            }

            geno_in_pheno <- length(intersect(rules_geno_fixed, rules_pheno_fixed)) / length(rules_geno_fixed)
            pheno_in_geno <- length(intersect(rules_pheno_fixed, rules_geno_fixed)) / length(rules_pheno_fixed)

            if (i == 1) {
                crs_geno <- geno_in_pheno
                crs_pheno <- pheno_in_geno
            } else {
                crs_geno <- (1/i) * ((i-1) * crs_geno + geno_in_pheno)
                crs_pheno <- (1/i) * ((i-1) * crs_pheno + pheno_in_geno)
            }

            crs_list_geno <- list.append(crs_list_geno, crs_geno)
            crs_list_pheno <- list.append(crs_list_pheno, crs_pheno)
            year_list <- list.append(year_list, year)
        }

        data.frame(
            year = unlist(year_list),
            geno_in_pheno = unlist(crs_list_geno),
            pheno_in_geno = unlist(crs_list_pheno),
            data_source = data_source
        )
    }

    crs_df1 <- process_dataset(genotype_df1, phenotype_df1, data_source1)
    crs_df2 <- process_dataset(genotype_df2, phenotype_df2, data_source2)

    crs_df <- rbind(crs_df1, crs_df2)

   plot1 <- ggplot(crs_df) +
    geom_line(aes(x = year, y = geno_in_pheno, color = paste(data_source, "Geno in Pheno")), linewidth = 2) +
    geom_point(aes(x = year, y = geno_in_pheno, color = paste(data_source, "Geno in Pheno")), size = 3) +
    geom_line(aes(x = year, y = pheno_in_geno, color = paste(data_source, "Pheno in Geno")), linewidth = 2) +
    geom_point(aes(x = year, y = pheno_in_geno, color = paste(data_source, "Pheno in Geno")), size = 3) +
    theme(
        text = element_text(size = 15),#a82c2c
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.8, "cm"),
        legend.box = "vertical",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    ) +
    labs(y = "Cumulative Proportion Captured", x = "Year", 
         title = "Cumulative Proportion of Rules Captured") +
    scale_x_continuous(breaks = round(seq(min(crs_df$year), max(crs_df$year), by = 2),1)) +
    scale_color_manual(values = c("blue", "red", "green", "purple"),
                       labels = c("Cecal Genotype in Phenotype", "Cecal Phenotype in Genotype", 
                                  "Retail Meats Genotype in Phenotype", "Retail Meats Phenotype in Genotype")) +
    ylim(0,1) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave(plot = plot1, filename = str_glue("all_combined/figures/shared_rules_analysis/between_indicators/crs/combined_{target}_captured_bidirectional_{rules_selected}.png"),
       width = 10, height = 8)
    return(plot1)
}


# Usage
genotype_df_retail <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_genotype.csv")
phenotype_df_retail <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_phenotype.csv")
genotype_df_cecal <- read.csv("cecal/cecal_wide_class_level_genotype.csv")
phenotype_df_cecal <- read.csv("cecal/cecal_wide_class_level_phenotype.csv")

 #Include only rows where the ID matches in both dataframes for each data source
phenotype_df_retail <- phenotype_df_retail[phenotype_df_retail$ID %in% genotype_df_retail$ID, ]
 genotype_df_retail <- genotype_df_retail[genotype_df_retail$ID %in% phenotype_df_retail$ID, ]

phenotype_df_cecal <- phenotype_df_cecal[phenotype_df_cecal$ID %in% genotype_df_cecal$ID, ]
genotype_df_cecal <- genotype_df_cecal[genotype_df_cecal$ID %in% phenotype_df_cecal$ID, ]

 gene_vs_phenotype_cumulative_percent_captured_multi(
    genotype_df1 = genotype_df_retail, 
    phenotype_df1 = phenotype_df_retail, 
    genotype_df2 = genotype_df_cecal, 
    phenotype_df2 = phenotype_df_cecal, 
    target = "rules", 
    rules_selected = "best", 
    cut_off = c(0.5, 0, 0.5, 0), 
    measures_used = c("cosine", "jaccard", "kulczynski", "support"), 
    data_source1 = "Retail_Meats",
    data_source2 = "Cecal"
)





################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
#####################################################   BETWEEN DATA SOURCES      ##############################################################################################
################################################################################################################################################################################    
################################################################################################################################################################################


# Function to calculate and plot cumulative percent captured between cecal and retail for both genotype and phenotype
cecal_vs_retail_cumulative_percent_captured_multi <- function(cecal_genotype_df, cecal_phenotype_df, retail_genotype_df, retail_phenotype_df, target, rules_selected = "best", cut_off, measures_used, data_source1, data_source2) {
    process_dataset_pair <- function(df1, df2, indicator, data_source1, data_source2) {
        crs_list_1in2 <- list()
        crs_list_2in1 <- list()
        year_list <- list()

        for (year in c(min(df1[, "Year"]):(max(df1[, "Year"])))) {
            i <- match(year, sort(unique(df1$Year)))
            if (rules_selected == "best") {
                rules1 <- implement_cuts(df = df1, resistance_indicator = indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
                rules2 <- implement_cuts(df = df2, resistance_indicator = indicator, target = target, cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
            } else if (rules_selected == "all") {
                rules1 <- get_rules_or_itemsets(df = df1, resistance_indicator = indicator, year = year, target = target, agg = FALSE)
                rules2 <- get_rules_or_itemsets(df = df2, resistance_indicator = indicator, year = year, target = target, agg = FALSE)
            }

            new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
            if (target == "rules") {
                rules1_fixed <- new("rules",
                    lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                    rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
                )
                rules2_fixed <- new("rules",
                    lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                    rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
                )
            } else {
                rules1_fixed <- arules::recode(rules1, itemLabels = new_itemLabels)
                rules2_fixed <- arules::recode(rules2, itemLabels = new_itemLabels)
            }

            d1in2 <- length(intersect(rules1_fixed, rules2_fixed)) / length(rules1_fixed)
            d2in1 <- length(intersect(rules2_fixed, rules1_fixed)) / length(rules2_fixed)

            if (i == 1) {
                crs_1in2 <- d1in2
                crs_2in1 <- d2in1
            } else {
                crs_1in2 <- (1/i) * ((i-1) * crs_1in2 + d1in2)
                crs_2in1 <- (1/i) * ((i-1) * crs_2in1 + d2in1)
            }

            crs_list_1in2 <- list.append(crs_list_1in2, crs_1in2)
            crs_list_2in1 <- list.append(crs_list_2in1, crs_2in1)
            year_list <- list.append(year_list, year)
        }

        data.frame(
            year = unlist(year_list),
            d1in2 = unlist(crs_list_1in2),
            d2in1 = unlist(crs_list_2in1),
            indicator = indicator
        )
    }

    crs_df_pheno <- process_dataset_pair(cecal_phenotype_df, retail_phenotype_df, "phenotype", data_source1, data_source2)
    crs_df_geno <- process_dataset_pair(cecal_genotype_df, retail_genotype_df, "genotype", data_source1, data_source2)

    crs_df <- rbind(crs_df_pheno, crs_df_geno)

    plot1 <- ggplot(crs_df) +
        geom_line(aes(x = year, y = d1in2, color = paste(indicator, data_source1, "in", data_source2)), linewidth = 2) +
        geom_point(aes(x = year, y = d1in2, color = paste(indicator, data_source1, "in", data_source2)), size = 3) +
        geom_line(aes(x = year, y = d2in1, color = paste(indicator, data_source2, "in", data_source1)), linewidth = 2) +
        geom_point(aes(x = year, y = d2in1, color = paste(indicator, data_source2, "in", data_source1)), size = 3) +
        theme(
            text = element_text(size = 15),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.8, "cm"),
            legend.box = "vertical",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
        ) +
        labs(y = "Cumulative Proportion Captured", x = "Year", 
             title = "Cumulative Proportion of Rules Captured Between Data Sources") +
        scale_x_continuous(breaks = round(seq(min(crs_df$year), max(crs_df$year), by = 2),1)) +
        scale_color_manual(values = c("blue", "red", "green", "purple"),
                           labels = c(paste("Genotype", data_source1, "in", data_source2), 
                                      paste("Genotype", data_source2, "in", data_source1),
                                      paste("Phenotype", data_source1, "in", data_source2), 
                                      paste("Phenotype", data_source2, "in", data_source1))) +
        ylim(0,1) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE))

    ggsave(plot = plot1, filename = str_glue("all_combined/figures/shared_rules_analysis/between_data_sources/crs/combined_{target}_captured_bidirectional_{rules_selected}.png"),
           width = 10, height = 8)
    return(plot1)
}




#only keep years in phenotype dfs that are present in both cecal and retail
# phenotype_df_cecal <- phenotype_df_cecal[phenotype_df_cecal$Year %in% phenotype_df_retail$Year, ]
# phenotype_df_retail <- phenotype_df_retail[phenotype_df_retail$Year %in% phenotype_df_cecal$Year, ]

#only keep years in genotype dfs that are present in both cecal and retail
# genotype_df_cecal <- genotype_df_cecal[genotype_df_cecal$Year %in% genotype_df_retail$Year, ]
# genotype_df_retail <- genotype_df_retail[genotype_df_retail$Year %in% genotype_df_cecal$Year, ]

# cecal_vs_retail_cumulative_percent_captured_multi(
    # cecal_genotype_df = genotype_df_cecal, 
    # cecal_phenotype_df = phenotype_df_cecal, 
    # retail_genotype_df = genotype_df_retail, 
    # retail_phenotype_df = phenotype_df_retail, 
    # target = "rules", 
    # rules_selected = "best", 
    # cut_off = c(0.5, 0, 0.5, 0), 
    # measures_used = c("cosine", "jaccard", "kulczynski", "support"), 
    # data_source1 = "Cecal",
    # data_source2 = "Retail Meats"
# )