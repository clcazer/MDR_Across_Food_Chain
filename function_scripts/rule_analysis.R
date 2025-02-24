

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

library(here)

source(here("function_scripts", "data_wrangling.R"))
source(here("function_scripts", "rule_mining_and_selection.R"))


# Function to create scatter plot comparing number of rules before and after implementing measure cut-offs
plot_all_vs_best_rules <- function(df, resistance_indicator, target, cut_off, measures_used, data_source, class_level = FALSE) {
    # Initialize lists to store data
    all_rules_length <- list()
    best_rules_length <- list()
    year_list <- list()
    sample_size <- list()

    # Get available years from the dataset
    available_years <- sort(unique(df[, "Year"]))
    
    # Iterate through each available year
    for (year in available_years) {
        tryCatch({
            # Get all rules for the current year
            all_rules <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator,
                                             year = year, target = target)
            
            # Only process if rules were found
            if (length(all_rules) > 0) {
                all_rules_length <- list.append(all_rules_length, length(all_rules))

                # Apply cut-offs to get best rules
                best_rules <- implement_cuts(df = df, resistance_indicator = resistance_indicator, 
                                          target = target, cut_off = cut_off, year = year, 
                                          measures_used = measures_used)
                
                best_rules_length <- list.append(best_rules_length, length(best_rules))
                year_list <- list.append(year_list, year)

                # Calculate sample size for the current year
                sample_size_df <- df[df$Year == year,]
                sample_size <- list.append(sample_size, length(unique(sample_size_df$ID)))
            }
        }, error = function(e) {
            warning(sprintf("Could not process year %d: %s", year, e$message))
        })
    }

    # Only create plot if we have data
    if (length(year_list) > 0) {
        # Create dataframe with collected data
        best_v_all_df <- data.frame(
            Sample_Size = unlist(sample_size), 
            Year = as.factor(unlist(year_list)), 
            All_Rules = unlist(all_rules_length), 
            Best_Rules = unlist(best_rules_length)
        )

        print(best_v_all_df)


        # save the df to have the actual numbers for reference
        if (class_level) {
        write.csv(best_v_all_df, file = str_glue("{data_source}/figures/pruning_results/{resistance_indicator}/{data_source}_{resistance_indicator}_Pruning_Results_CLASS.csv"), row.names = FALSE)
        } else {
        write.csv(best_v_all_df, file = str_glue("{data_source}/figures/pruning_results/{resistance_indicator}/{data_source}_{resistance_indicator}_Pruning_Results.csv"), row.names = FALSE)
        }


        # Create and save the plot
        plot <- ggplot(best_v_all_df, aes(x = All_Rules, y = Best_Rules, col = Year)) + 
            geom_point(size = 7) +
            theme(text = element_text(size = 20), 
                  plot.title = element_text(hjust = 0.5),
                  panel.background = element_rect(fill = 'white', colour = 'black')) +
            labs(y = str_glue("Best Rules"), 
                 x = str_glue("All Rules"), 
                 title = str_glue("Rule Pruning Results"))

        # Save the plot with appropriate filename
        filename <- if (class_level) {
            str_glue("{data_source}/figures/pruning_results/{resistance_indicator}/{data_source}_{resistance_indicator}_{target}_Rule_Pruning_Results_CLASS.png")
        } else {
            str_glue("{data_source}/figures/pruning_results/{resistance_indicator}/{data_source}_{resistance_indicator}_{target}_Rule_Pruning_Results.png")
        }
        
        ggsave(filename)
        return(plot)
    } else {
        warning("No valid data found for plotting")
        return(NULL)
    }
}








################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
#####################################################   WITHIN INDICATORS      ################################################################################################
################################################################################################################################################################################    
################################################################################################################################################################################


# Function to calculate and plot cumulative proportion captured for both genotype and phenotype across two data sources
cumulative_proportion_captured_multi <- function(datasets_geno, datasets_pheno, class_level = FALSE, target, rules_selected = "best", cut_off, measures_used, data_partition) {
    all_crs_df <- data.frame()

    # Process both genotype and phenotype datasets
    for (resistance_indicator in c("genotype", "phenotype")) {
        current_datasets <- if(resistance_indicator == "genotype") datasets_geno else datasets_pheno
        
        for (dataset_name in names(current_datasets)) {
            df <- current_datasets[[dataset_name]]
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

                # Fix item encodings
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

            crs_df <- data.frame(
                crs = unlist(crs_list), 
                year = unlist(year_list), 
                dataset = dataset_name,
                indicator = resistance_indicator
            )
            all_crs_df <- rbind(all_crs_df, crs_df)
        }
    }

    
    class_prefix <- if(class_level) "Class " else ""
    title_prefix <- if(class_level) "Class Level " else "Individual Antimicrobial/Gene "

    # Create a mapping of dataset names to their visual properties
    if (!class_level) {
        visual_mapping <- list(
            "Cecal Genotype" = list(color = "#0077BB", linetype = "solid", shape = 16, label = "Cecal Genotype"),
            "Retail Meats Genotype" = list(color = "#009E73", linetype = "solid", shape = 17, label = "Retail Meats Genotype"),
            "NAHLN Genotype" = list(color = "#EE7733", linetype = "solid", shape = 15, label = "NAHLN Genotype"),
            "Cecal Phenotype" = list(color = "#0077BB", linetype = "dashed", shape = 16, label = "Cecal Phenotype"),
            "Retail Meats Phenotype" = list(color = "#009E73", linetype = "dashed", shape = 17, label = "Retail Meats Phenotype"),
            "NAHLN Phenotype" = list(color = "#EE7733", linetype = "dashed", shape = 15, label = "NAHLN Phenotype")
        )
    } else {
        visual_mapping <- list(
            "Cecal Class Genotype" = list(color = "#0077BB", linetype = "solid", shape = 16, label = "Cecal Class Genotype"),
            "Retail Meats Class Genotype" = list(color = "#009E73", linetype = "solid", shape = 17, label = "Retail Meats Class Genotype"),
            "NAHLN Class Genotype" = list(color = "#EE7733", linetype = "solid", shape = 15, label = "NAHLN Class Genotype"),
            "Cecal Class Phenotype" = list(color = "#0077BB", linetype = "dashed", shape = 16, label = "Cecal Class Phenotype"),
            "Retail Meats Class Phenotype" = list(color = "#009E73", linetype = "dashed", shape = 17, label = "Retail Meats Class Phenotype"),
            "NAHLN Class Phenotype" = list(color = "#EE7733", linetype = "dashed", shape = 15, label = "NAHLN Class Phenotype")
        )
    }

    # Create vectors for scale_*_manual
    color_values <- sapply(visual_mapping, function(x) x$color)
    linetype_values <- sapply(visual_mapping, function(x) x$linetype)
    shape_values <- sapply(visual_mapping, function(x) x$shape)
    labels <- sapply(visual_mapping, function(x) x$label)

    # Create the plot using the mapping
    plot1 <- ggplot(all_crs_df, aes(x = year, y = crs, color = dataset, linetype = dataset, shape = dataset)) +
        geom_line(linewidth = 1.5) +
        geom_point(size = 4) +
        theme(
            text = element_text(size = 25),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 18),
            legend.key.size = unit(1.2, "cm"),
            legend.box = "vertical",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
        ) +
        labs(y = "Cumulative Proportion Captured", x = "Year", 
             title = str_glue("{title_prefix}Cumulative Rule Set Comparison")) +
        scale_x_continuous(breaks = round(seq(min(all_crs_df$year), max(all_crs_df$year), by = 2),1)) +
        scale_color_manual(
            values = color_values,
            labels = labels,
            breaks = names(visual_mapping)
        ) +
        scale_linetype_manual(
            values = linetype_values,
            labels = labels,
            breaks = names(visual_mapping)
        ) +
        scale_shape_manual(
            values = shape_values,
            labels = labels,
            breaks = names(visual_mapping)
        ) +
        ylim(0,1) +
        guides(
            color = guide_legend(nrow = 3, byrow = TRUE),
            linetype = guide_legend(nrow = 3, byrow = TRUE),
            shape = guide_legend(nrow = 3, byrow = TRUE)
        )

    filename_suffix <- if(class_level) "_CLASS" else ""
    ggsave(plot = plot1, 
           filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/within_indicators/multi_combined_{target}_c_prop_captured_{rules_selected}{filename_suffix}.png"),
           width = 10, height = 8)
    
    return(plot1)
}




cumulative_proportion_captured_multi_v2 <- function(datasets_geno, datasets_pheno, class_level = FALSE, target, rules_selected = "best", cut_off, measures_used, data_partition) {
    all_crs_df <- data.frame()

    # Process both genotype and phenotype datasets
    for (resistance_indicator in c("genotype", "phenotype")) {
        current_datasets <- if(resistance_indicator == "genotype") datasets_geno else datasets_pheno
        
        print(paste("Processing", resistance_indicator, "datasets"))
        print(paste("Number of datasets:", length(current_datasets)))
        
        for (dataset_name in names(current_datasets)) {
            print(paste("Processing dataset:", dataset_name))
            df <- current_datasets[[dataset_name]]
            crs_list <- list()
            year_list <- list()
            
            # Get sorted unique years for this dataset
            available_years <- sort(unique(df$Year))
            print(paste("Available years:", paste(available_years, collapse = ", ")))
            
            # Process each year except the last one
            for (i in 1:(length(available_years)-1)) {
                current_year <- available_years[i]
                next_year <- available_years[i + 1]
                
                print(paste("Processing year:", current_year, "with next year:", next_year))
                
                tryCatch({
                    if (rules_selected == "best") {
                        rules1 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, 
                                              target = target, cut_off = cut_off, 
                                              measures_used = measures_used, year = current_year)
                        rules2 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, 
                                              target = target, cut_off = cut_off, 
                                              measures_used = measures_used, year = next_year)
                    } else if (rules_selected == "all") {
                        rules1 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, 
                                                      year = current_year, target = target)
                        rules2 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, 
                                                      year = next_year, target = target)
                    }
                    
                    print(paste("Number of rules - Year", current_year, ":", length(rules1)))
                    print(paste("Number of rules - Year", next_year, ":", length(rules2)))
                    
                    # Only proceed if both rule sets have items
                    if (length(rules1) > 0 && length(rules2) > 0) {
                        new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
                        
                        if (target == "rules") {
                            # Check if both rules have LHS and RHS before recoding
                            if (length(lhs(rules1)) > 0 && length(rhs(rules1)) > 0 && 
                                length(lhs(rules2)) > 0 && length(rhs(rules2)) > 0) {
                                
                                rules1_fixed <- new("rules",
                                    lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                                    rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
                                )
                                rules2_fixed <- new("rules",
                                    lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                                    rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
                                )
                                
                                intersection_size <- length(intersect(rules1_fixed, rules2_fixed))
                                if (i == 1) {
                                    crs <- intersection_size/length(rules1_fixed)
                                } else {
                                    crs <- (1/i)*((i-1)*crs + intersection_size/length(rules1_fixed))
                                }
                                
                                crs_list <- list.append(crs_list, crs)
                                year_list <- list.append(year_list, current_year)
                            }
                        } else {
                            # Handle itemsets case
                            rules1_fixed <- arules::recode(rules1, itemLabels = new_itemLabels)
                            rules2_fixed <- arules::recode(rules2, itemLabels = new_itemLabels)
                            
                            intersection_size <- length(intersect(rules1_fixed, rules2_fixed))
                            if (i == 1) {
                                crs <- intersection_size/length(rules1_fixed)
                            } else {
                                crs <- (1/i)*((i-1)*crs + intersection_size/length(rules1_fixed))
                            }
                            
                            crs_list <- list.append(crs_list, crs)
                            year_list <- list.append(year_list, current_year)
                        }
                    }
                }, error = function(e) {
                    print(paste("Error processing year", current_year, ":", e$message))
                })
            }

            # Only create dataframe if we have results
            if (length(crs_list) > 0) {
                crs_df <- data.frame(
                    crs = unlist(crs_list), 
                    year = unlist(year_list), 
                    dataset = dataset_name,
                    indicator = resistance_indicator
                )
                all_crs_df <- rbind(all_crs_df, crs_df)
            }
        }
    }

    # Only proceed with plotting if we have data
    if (nrow(all_crs_df) > 0) {
        class_prefix <- if(class_level) "Class " else ""
        title_prefix <- if(class_level) "Class Level " else "Individual Antimicrobial/Gene "

        # Create a mapping of dataset names to their visual properties
        if (!class_level) {
            visual_mapping <- list(
                "Cecal Genotype" = list(color = "#0077BB", linetype = "solid", shape = 16, label = "Cecal Genotype"),
                "Retail Meats Genotype" = list(color = "#009E73", linetype = "solid", shape = 17, label = "Retail Meats Genotype"),
                "NAHLN Genotype" = list(color = "#EE7733", linetype = "solid", shape = 15, label = "NAHLN Genotype"),
                "Cecal Phenotype" = list(color = "#0077BB", linetype = "dashed", shape = 16, label = "Cecal Phenotype"),
                "Retail Meats Phenotype" = list(color = "#009E73", linetype = "dashed", shape = 17, label = "Retail Meats Phenotype"),
                "NAHLN Phenotype" = list(color = "#EE7733", linetype = "dashed", shape = 15, label = "NAHLN Phenotype")
            )
        } else {
            visual_mapping <- list(
                "Cecal Class Genotype" = list(color = "#0077BB", linetype = "solid", shape = 16, label = "Cecal Class Genotype"),
                "Retail Meats Class Genotype" = list(color = "#009E73", linetype = "solid", shape = 17, label = "Retail Meats Class Genotype"),
                "NAHLN Class Genotype" = list(color = "#EE7733", linetype = "solid", shape = 15, label = "NAHLN Class Genotype"),
                "Cecal Class Phenotype" = list(color = "#0077BB", linetype = "dashed", shape = 16, label = "Cecal Class Phenotype"),
                "Retail Meats Class Phenotype" = list(color = "#009E73", linetype = "dashed", shape = 17, label = "Retail Meats Class Phenotype"),
                "NAHLN Class Phenotype" = list(color = "#EE7733", linetype = "dashed", shape = 15, label = "NAHLN Class Phenotype")
            )
        }

        # Create vectors for scale_*_manual
        color_values <- sapply(visual_mapping, function(x) x$color)
        linetype_values <- sapply(visual_mapping, function(x) x$linetype)
        shape_values <- sapply(visual_mapping, function(x) x$shape)
        labels <- sapply(visual_mapping, function(x) x$label)

        # Create the plot using the mapping
        plot1 <- ggplot(all_crs_df, aes(x = year, y = crs, color = dataset, linetype = dataset, shape = dataset)) +
            geom_line(linewidth = 1.5) +
            geom_point(size = 4) +
            theme(
                text = element_text(size = 25),
                plot.title = element_text(hjust = 0.5),
                panel.background = element_rect(fill = 'white', colour = 'black'),
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.text = element_text(size = 18),
                legend.key.size = unit(1.2, "cm"),
                legend.box = "vertical",
                legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
            ) +
            labs(y = "Cumulative Proportion Captured", x = "Year", 
                 title = str_glue("{title_prefix}Cumulative Rule Set Comparison")) +
            scale_x_continuous(breaks = round(seq(min(all_crs_df$year), max(all_crs_df$year), by = 2),1)) +
            scale_color_manual(
                values = color_values,
                labels = labels,
                breaks = names(visual_mapping)
            ) +
            scale_linetype_manual(
                values = linetype_values,
                labels = labels,
                breaks = names(visual_mapping)
            ) +
            scale_shape_manual(
                values = shape_values,
                labels = labels,
                breaks = names(visual_mapping)
            ) +
            ylim(0,1) +
            guides(
                color = guide_legend(nrow = 3, byrow = TRUE),
                linetype = guide_legend(nrow = 3, byrow = TRUE),
                shape = guide_legend(nrow = 3, byrow = TRUE)
            )

        filename_suffix <- if(class_level) "_CLASS" else ""
        ggsave(plot = plot1, 
               filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/within_indicators/multi_combined_{target}_c_prop_captured_{rules_selected}{filename_suffix}.png"),
               width = 10, height = 8)
        
        return(plot1)
    } else {
        stop("No data available for plotting")
    }
}


#########NON CUMULATIVE WITHIN INDICATORS COMPARISONS
proportion_captured_multi_v2 <- function(datasets_geno, datasets_pheno, class_level = FALSE, target, rules_selected = "best", cut_off, measures_used, data_partition) {
    all_prop_df <- data.frame()

    # Process both genotype and phenotype datasets
    for (resistance_indicator in c("genotype", "phenotype")) {
        current_datasets <- if(resistance_indicator == "genotype") datasets_geno else datasets_pheno
        
        print(paste("Processing", resistance_indicator, "datasets"))
        print(paste("Number of datasets:", length(current_datasets)))
        
        for (dataset_name in names(current_datasets)) {
            print(paste("Processing dataset:", dataset_name))
            df <- current_datasets[[dataset_name]]
            prop_list <- list()
            year_list <- list()
            
            # Get sorted unique years for this dataset
            available_years <- sort(unique(df$Year))
            print(paste("Available years:", paste(available_years, collapse = ", ")))
            
            # Process each year except the last one
            for (i in 1:(length(available_years)-1)) {
                current_year <- available_years[i]
                next_year <- available_years[i + 1]
                
                print(paste("Processing year:", current_year, "with next year:", next_year))
                
                tryCatch({
                    if (rules_selected == "best") {
                        rules1 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, 
                                              target = target, cut_off = cut_off, 
                                              measures_used = measures_used, year = current_year)
                        rules2 <- implement_cuts(df = df, resistance_indicator = resistance_indicator, 
                                              target = target, cut_off = cut_off, 
                                              measures_used = measures_used, year = next_year)
                    } else if (rules_selected == "all") {
                        rules1 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, 
                                                      year = current_year, target = target)
                        rules2 <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, 
                                                      year = next_year, target = target)
                    }
                    
                    print(paste("Number of rules - Year", current_year, ":", length(rules1)))
                    print(paste("Number of rules - Year", next_year, ":", length(rules2)))
                    
                    # Only proceed if both rule sets have items
                    if (length(rules1) > 0 && length(rules2) > 0) {
                        new_itemLabels <- union(itemLabels(rules1), itemLabels(rules2))
                        
                        if (target == "rules") {
                            # Check if both rules have LHS and RHS before recoding
                            if (length(lhs(rules1)) > 0 && length(rhs(rules1)) > 0 && 
                                length(lhs(rules2)) > 0 && length(rhs(rules2)) > 0) {
                                
                                rules1_fixed <- new("rules",
                                    lhs = arules::recode(lhs(rules1), itemLabels = new_itemLabels), 
                                    rhs = arules::recode(rhs(rules1), itemLabels = new_itemLabels)
                                )
                                rules2_fixed <- new("rules",
                                    lhs = arules::recode(lhs(rules2), itemLabels = new_itemLabels), 
                                    rhs = arules::recode(rhs(rules2), itemLabels = new_itemLabels)
                                )
                                
                                # Calculate proportion directly without cumulative averaging
                                prop <- length(intersect(rules1_fixed, rules2_fixed))/length(rules1_fixed)
                                
                                prop_list <- list.append(prop_list, prop)
                                year_list <- list.append(year_list, current_year)
                            }
                        } else {
                            # Handle itemsets case
                            rules1_fixed <- arules::recode(rules1, itemLabels = new_itemLabels)
                            rules2_fixed <- arules::recode(rules2, itemLabels = new_itemLabels)
                            
                            # Calculate proportion directly without cumulative averaging
                            prop <- length(intersect(rules1_fixed, rules2_fixed))/length(rules1_fixed)
                            
                            prop_list <- list.append(prop_list, prop)
                            year_list <- list.append(year_list, current_year)
                        }
                    }
                }, error = function(e) {
                    print(paste("Error processing year", current_year, ":", e$message))
                })
            }

            # Only create dataframe if we have results
            if (length(prop_list) > 0) {
                prop_df <- data.frame(
                    proportion = unlist(prop_list), 
                    year = unlist(year_list), 
                    dataset = dataset_name,
                    indicator = resistance_indicator
                )
                all_prop_df <- rbind(all_prop_df, prop_df)
            }
        }
    }

    # Only proceed with plotting if we have data
    if (nrow(all_prop_df) > 0) {
        class_prefix <- if(class_level) "Class " else ""
        title_prefix <- if(class_level) "Class Level " else "Individual Antimicrobial/Gene "

        # Create a mapping of dataset names to their visual properties
        if (!class_level) {
            visual_mapping <- list(
                "Cecal Genotype" = list(color = "#0077BB", linetype = "solid", shape = 16, label = "Cecal Genotype"),
                "Retail Meats Genotype" = list(color = "#009E73", linetype = "solid", shape = 17, label = "Retail Meats Genotype"),
                "NAHLN Genotype" = list(color = "#EE7733", linetype = "solid", shape = 15, label = "NAHLN Genotype"),
                "Cecal Phenotype" = list(color = "#0077BB", linetype = "dashed", shape = 16, label = "Cecal Phenotype"),
                "Retail Meats Phenotype" = list(color = "#009E73", linetype = "dashed", shape = 17, label = "Retail Meats Phenotype"),
                "NAHLN Phenotype" = list(color = "#EE7733", linetype = "dashed", shape = 15, label = "NAHLN Phenotype")
            )
        } else {
            visual_mapping <- list(
                "Cecal Class Genotype" = list(color = "#0077BB", linetype = "solid", shape = 16, label = "Cecal Class Genotype"),
                "Retail Meats Class Genotype" = list(color = "#009E73", linetype = "solid", shape = 17, label = "Retail Meats Class Genotype"),
                "NAHLN Class Genotype" = list(color = "#EE7733", linetype = "solid", shape = 15, label = "NAHLN Class Genotype"),
                "Cecal Class Phenotype" = list(color = "#0077BB", linetype = "dashed", shape = 16, label = "Cecal Class Phenotype"),
                "Retail Meats Class Phenotype" = list(color = "#009E73", linetype = "dashed", shape = 17, label = "Retail Meats Class Phenotype"),
                "NAHLN Class Phenotype" = list(color = "#EE7733", linetype = "dashed", shape = 15, label = "NAHLN Class Phenotype")
            )
        }

        # Create vectors for scale_*_manual
        color_values <- sapply(visual_mapping, function(x) x$color)
        linetype_values <- sapply(visual_mapping, function(x) x$linetype)
        shape_values <- sapply(visual_mapping, function(x) x$shape)
        labels <- sapply(visual_mapping, function(x) x$label)

        # Create the plot using the mapping
        plot1 <- ggplot(all_prop_df, aes(x = year, y = proportion, color = dataset, linetype = dataset, shape = dataset)) +
            geom_line(linewidth = 1.5) +
            geom_point(size = 4) +
            theme(
                text = element_text(size = 25),
                plot.title = element_text(hjust = 0.5),
                panel.background = element_rect(fill = 'white', colour = 'black'),
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.text = element_text(size = 18),
                legend.key.size = unit(1.2, "cm"),
                legend.box = "vertical",
                legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
            ) +
            labs(y = "Proportion Captured", x = "Year", 
                 title = str_glue("{title_prefix}Rule Set Comparison")) +
            scale_x_continuous(breaks = round(seq(min(all_prop_df$year), max(all_prop_df$year), by = 2),1)) +
            scale_color_manual(
                values = color_values,
                labels = labels,
                breaks = names(visual_mapping)
            ) +
            scale_linetype_manual(
                values = linetype_values,
                labels = labels,
                breaks = names(visual_mapping)
            ) +
            scale_shape_manual(
                values = shape_values,
                labels = labels,
                breaks = names(visual_mapping)
            ) +
            ylim(0,1) +
            guides(
                color = guide_legend(nrow = 3, byrow = TRUE),
                linetype = guide_legend(nrow = 3, byrow = TRUE),
                shape = guide_legend(nrow = 3, byrow = TRUE)
            )

        filename_suffix <- if(class_level) "_CLASS" else ""
        ggsave(plot = plot1, 
               filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/within_indicators/multi_combined_{target}_prop_captured_{rules_selected}{filename_suffix}.png"),
               width = 10, height = 8)
        
        return(plot1)
    } else {
        stop("No data available for plotting")
    }
}



# #load in the data

# #cecal and retail meats
# cecal_phenotype_df <- read.csv("cecal/NEW_cecal_wide_resStatus_phenotype.csv")
# cecal_Class_phenotype_df <- read.csv("cecal/NEW_cecal_wide_class_level_phenotype.csv")

# cecal_Family_genotype_df <- read.csv("cecal/cecal_wide_corrected_genotype.csv")
# cecal_Class_genotype_df <- read.csv("cecal/cecal_wide_class_level_genotype.csv")

# Retail_Meats_phenotype_df <- read.csv("Retail_Meats/NEW_Retail_Meats_wide_resStatus_phenotype.csv")
# Retail_Meats_Class_phenotype_df <- read.csv("Retail_Meats/NEW_Retail_Meats_wide_class_level_phenotype.csv")

# Retail_Meats_Family_genotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_corrected_genotype.csv")
# Retail_Meats_Class_genotype_df <- read.csv("Retail_Meats/Retail_Meats_wide_class_level_genotype.csv")

# #cecal and retail meats 2017 onward
# cecal_2017_phenotype_df <- read.csv("cecal_2017Onward/NEW_cecal_2017Onward_wide_resStatus_phenotype.csv")
# cecal_2017_Class_phenotype_df <- read.csv("cecal_2017Onward/NEW_cecal_2017Onward_wide_class_level_phenotype.csv")

# cecal_2017_Family_genotype_df <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_corrected_genotype.csv")
# cecal_2017_Class_genotype_df <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_class_level_genotype.csv")

# Retail_2017_Meats_phenotype_df <- read.csv("Retail_Meats_2017Onward/NEW_Retail_Meats_2017Onward_wide_resStatus_phenotype.csv")
# Retail_2017_Meats_Class_phenotype_df <- read.csv("Retail_Meats_2017Onward/NEW_Retail_Meats_2017Onward_wide_class_level_phenotype.csv")

# Retail_2017_Meats_Family_genotype_df <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_corrected_genotype.csv")
# Retail_2017_Meats_Class_genotype_df <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_class_level_genotype.csv")

# #NAHLN
# NAHLN_phenotype_df <- read.csv("NAHLN/NAHLN_wide_resStatus_phenotype.csv")
# NAHLN_Class_phenotype_df <- read.csv("NAHLN/NAHLN_wide_class_level_phenotype.csv")

# NAHLN_Family_genotype_df <- read.csv("NAHLN/NAHLN_wide_corrected_genotype.csv")
# NAHLN_Class_genotype_df <-read.csv("NAHLN/NAHLN_wide_class_level_genotype.csv")


# phenotype_datasets <- list(
#     "Cecal Phenotype" = cecal_phenotype_df,
#     "Retail Meats Phenotype" = Retail_Meats_phenotype_df,
#     "NAHLN Phenotype" = NAHLN_phenotype_df
# )

# class_phenotype_datasets <- list(
#     "Cecal Class Phenotype" = cecal_Class_phenotype_df,
#     "Retail Meats Class Phenotype" = Retail_Meats_Class_phenotype_df,
#     "NAHLN Class Phenotype" = NAHLN_Class_phenotype_df
# )

# genotype_datasets <- list(
#     "Cecal Genotype" = cecal_Family_genotype_df,
#     "Retail Meats Genotype" = Retail_Meats_Family_genotype_df,
#     "NAHLN Genotype" = NAHLN_Family_genotype_df
# )

# class_genotype_datasets <- list(
#     "Cecal Class Genotype" = cecal_Class_genotype_df,
#     "Retail Meats Class Genotype" = Retail_Meats_Class_genotype_df,
#     "NAHLN Class Genotype" = NAHLN_Class_genotype_df
# )

# data_partition = "full_date_range"

#  nonclass_plot <- proportion_captured_multi_v2(
#     datasets_geno = genotype_datasets,
#     datasets_pheno = phenotype_datasets,
#     class_level = FALSE,
#     target = "rules",
#     rules_selected = "best",
#     cut_off = c(0.5, 0, 0.5, 0),
#     measures_used = c("cosine", "jaccard", "kulczynski", "support"),
#     data_partition = data_partition
# )

#  class_plot <- proportion_captured_multi_v2(
#     datasets_geno = class_genotype_datasets,
#     datasets_pheno = class_phenotype_datasets,
#     class_level = TRUE,
#     target = "rules",
#     rules_selected = "best",
#     cut_off = c(0.5, 0, 0.5, 0),
#     measures_used = c("cosine", "jaccard", "kulczynski", "support"),
#     data_partition = data_partition
# )


# print("DONE")




################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
#####################################################   BETWEEN INDICATORS      ################################################################################################
################################################################################################################################################################################    
################################################################################################################################################################################




# Function to calculate and plot cumulative percent captured between genotype and phenotype for two data sources
gene_vs_phenotype_cumulative_percent_captured_multi <- function(genotype_df1, phenotype_df1, genotype_df2, phenotype_df2, target, rules_selected = "best", cut_off, measures_used, data_source1, data_source2, data_partition, title) {
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

    # Reshape the data for consistent aesthetics
   crs_df_long <- data.frame(
       year = rep(crs_df$year, 2),
       value = c(crs_df$geno_in_pheno, crs_df$pheno_in_geno),
       data_source = rep(crs_df$data_source, 2),
       type = rep(c("Geno in Pheno", "Pheno in Geno"), each = nrow(crs_df))
   )
   crs_df_long$interaction <- paste(crs_df_long$data_source, crs_df_long$type)
    # Create visual mapping
   visual_mapping <- list(
       "Cecal Geno in Pheno" = list(color = "#0077BB", linetype = "solid", shape = 16, label = "Cecal Genotype in Phenotype"),
       "Cecal Pheno in Geno" = list(color = "#0077BB", linetype = "dashed", shape = 16, label = "Cecal Phenotype in Genotype"),
       "Retail Meats Geno in Pheno" = list(color = "#009E73", linetype = "solid", shape = 17, label = "Retail Meats Genotype in Phenotype"),
       "Retail Meats Pheno in Geno" = list(color = "#009E73", linetype = "dashed", shape = 17, label = "Retail Meats Phenotype in Genotype")
   )
    # Create vectors for scale_*_manual
   color_values <- sapply(visual_mapping, function(x) x$color)
   linetype_values <- sapply(visual_mapping, function(x) x$linetype)
   shape_values <- sapply(visual_mapping, function(x) x$shape)
   labels <- sapply(visual_mapping, function(x) x$label)
    plot1 <- ggplot(crs_df_long, aes(x = year, y = value, color = interaction, linetype = interaction, shape = interaction)) +
       geom_line(linewidth = 1.5) +
       geom_point(size = 4) +
       theme(
           text = element_text(size = 25),
           plot.title = element_text(hjust = 0.5),
           panel.background = element_rect(fill = 'white', colour = 'black'),
           legend.position = "bottom",
           legend.title = element_blank(),
           legend.text = element_text(size = 18),
           legend.key.size = unit(1.2, "cm"),
           legend.box = "vertical",
           legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
       ) +
       labs(y = "Cumulative Proportion Captured", x = "Year", 
            title = str_glue({title})) + 
       scale_x_continuous(breaks = round(seq(min(crs_df$year), max(crs_df$year), by = 2),1)) +
       scale_color_manual(
           values = color_values,
           labels = labels,
           breaks = names(visual_mapping)
       ) +
       scale_linetype_manual(
           values = linetype_values,
           labels = labels,
           breaks = names(visual_mapping)
       ) +
       scale_shape_manual(
           values = shape_values,
           labels = labels,
           breaks = names(visual_mapping)
       ) +
       ylim(0,1) +
       guides(
           color = guide_legend(nrow = 2, byrow = TRUE),
           linetype = guide_legend(nrow = 2, byrow = TRUE),
           shape = guide_legend(nrow = 2, byrow = TRUE)
       )

ggsave(plot = plot1, filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/between_indicators/c_proportion_captured/combined_{target}_captured_bidirectional_{rules_selected}.png"),
       width = 10, height = 8)
    return(plot1)
}






# Function to calculate and plot cumulative percent captured between genotype and phenotype for three data sources
gene_vs_phenotype_cumulative_percent_captured_multi_v2 <- function(datasets_geno, datasets_pheno, target, rules_selected = "best", 
                                                                  cut_off, measures_used, data_partition, title) {
    # Debug print
    print("Starting function")
    print(paste("Number of genotype datasets:", length(datasets_geno)))
    print(paste("Number of phenotype datasets:", length(datasets_pheno)))
    print("Dataset names:")
    print(names(datasets_geno))
    print(names(datasets_pheno))

    process_dataset <- function(genotype_df, phenotype_df, data_source) {
        print(paste("Processing dataset for:", data_source))
        crs_list_geno <- list()
        crs_list_pheno <- list()
        year_list <- list()

        # Get common years between genotype and phenotype datasets
        geno_years <- sort(unique(genotype_df$Year))
        pheno_years <- sort(unique(phenotype_df$Year))
        available_years <- sort(intersect(geno_years, pheno_years))
        
        print(paste("Available years:", paste(available_years, collapse = ", ")))

        for (year in available_years) {
            print(paste("Processing year:", year))
            i <- match(year, available_years)
            
            tryCatch({
                # Get rules
                if (rules_selected == "best") {
                    rules_pheno <- implement_cuts(df = phenotype_df, resistance_indicator = "phenotype", 
                                               target = target, cut_off = cut_off, 
                                               measures_used = measures_used, year = year, agg = FALSE)
                    rules_geno <- implement_cuts(df = genotype_df, resistance_indicator = "genotype", 
                                              target = target, cut_off = cut_off, 
                                              measures_used = measures_used, year = year, agg = FALSE)
                } else if (rules_selected == "all") {
                    rules_pheno <- get_rules_or_itemsets(df = phenotype_df, resistance_indicator = "phenotype", 
                                                       year = year, target = target, agg = FALSE)
                    rules_geno <- get_rules_or_itemsets(df = genotype_df, resistance_indicator = "genotype", 
                                                      year = year, target = target, agg = FALSE)
                }

                print(paste("Number of rules - Pheno:", length(rules_pheno), "Geno:", length(rules_geno)))

                # Only process if both rule sets have items
                if (length(rules_pheno) > 0 && length(rules_geno) > 0) {
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
                } else {
                    print(paste("Skipping year", year, "- insufficient rules"))
                }
            }, error = function(e) {
                print(paste("Error processing year", year, ":", e$message))
            })
        }

        # Only create dataframe if we have results
        if (length(year_list) > 0) {
            print("Creating data frame for dataset")
            result_df <- data.frame(
                year = c(unlist(year_list), unlist(year_list)),
                value = c(unlist(crs_list_geno), unlist(crs_list_pheno)),
                type = rep(c("Geno in Pheno", "Pheno in Geno"), each = length(year_list)),
                data_source = data_source
            )
            print(paste("Number of rows in result:", nrow(result_df)))
            return(result_df)
        } else {
            print("No valid data found for this dataset pair")
            return(NULL)
        }
    }

    all_crs_df <- data.frame()
    
    # Process each dataset pair
    for (dataset_name in names(datasets_geno)) {
        print(paste("Processing dataset pair:", dataset_name))
        base_name <- gsub(" Genotype", "", dataset_name)
        pheno_name <- paste(base_name, "Phenotype")
        
        crs_df <- process_dataset(
            genotype_df = datasets_geno[[dataset_name]],
            phenotype_df = datasets_pheno[[pheno_name]],
            data_source = base_name
        )
        
        all_crs_df <- rbind(all_crs_df, crs_df)
    }
    
    # Create the interaction variable for plotting
    all_crs_df$interaction <- paste(all_crs_df$data_source, all_crs_df$type)

 # Create visual mapping
   visual_mapping <- list(
       "Cecal Geno in Pheno" = list(color = "#0077BB", linetype = "solid", shape = 16, label = "Cecal Genotype in Phenotype"),
       "Cecal Pheno in Geno" = list(color = "#0077BB", linetype = "dashed", shape = 16, label = "Cecal Phenotype in Genotype"),
       "Retail Meats Geno in Pheno" = list(color = "#009E73", linetype = "solid", shape = 17, label = "Retail Meats Genotype in Phenotype"),
       "Retail Meats Pheno in Geno" = list(color = "#009E73", linetype = "dashed", shape = 17, label = "Retail Meats Phenotype in Genotype"),
       "NAHLN Geno in Pheno" = list(color = "#EE7733", linetype = "solid", shape = 15, label = "NAHLN Genotype in Phenotype"),
       "NAHLN Pheno in Geno" = list(color = "#EE7733", linetype = "dashed", shape = 15, label = "NAHLN Phenotype in Genotype")
   )
    # Create vectors for scale_*_manual
   color_values <- sapply(visual_mapping, function(x) x$color)
   linetype_values <- sapply(visual_mapping, function(x) x$linetype)
   shape_values <- sapply(visual_mapping, function(x) x$shape)
   labels <- sapply(visual_mapping, function(x) x$label)
    plot1 <- ggplot(all_crs_df, aes(x = year, y = value, color = interaction, linetype = interaction, shape = interaction)) +
       geom_line(linewidth = 1.5) +
       geom_point(size = 4) +
       theme(
           text = element_text(size = 25),
           plot.title = element_text(hjust = 0.5),
           panel.background = element_rect(fill = 'white', colour = 'black'),
           legend.position = "bottom",
           legend.title = element_blank(),
           legend.text = element_text(size = 15),
           legend.key.size = unit(1.2, "cm"),
           legend.box = "vertical",
           legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
       ) +
       labs(y = "Cumulative Proportion Captured", x = "Year", 
            title = str_glue({title})) + 
       scale_x_continuous(breaks = round(seq(min(all_crs_df$year), max(all_crs_df$year), by = 2),1)) +
       scale_color_manual(
           values = color_values,
           labels = labels,
           breaks = names(visual_mapping)
       ) +
       scale_linetype_manual(
           values = linetype_values,
           labels = labels,
           breaks = names(visual_mapping)
       ) +
       scale_shape_manual(
           values = shape_values,
           labels = labels,
           breaks = names(visual_mapping)
       ) +
       ylim(0,1) +
       guides(
           color = guide_legend(nrow = 3, byrow = TRUE),
           linetype = guide_legend(nrow = 3, byrow = TRUE),
           shape = guide_legend(nrow = 3, byrow = TRUE)
       )

    ggsave(plot = plot1, 
           filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/between_indicators/c_proportion_captured/combined_{target}_captured_bidirectional_{rules_selected}.png"),
           width = 10, height = 8)
    
    return(plot1)
}






# Function to calculate and plot cumulative percent captured between genotype and phenotype for three data sources
gene_vs_phenotype_cumulative_percent_captured_multi_v3 <- function(datasets_geno, datasets_pheno, target, rules_selected = "best", 
                                                                  cut_off, measures_used, data_partition, title) {
    # Debug print
    print("Starting function")
    print(paste("Number of genotype datasets:", length(datasets_geno)))
    print(paste("Number of phenotype datasets:", length(datasets_pheno)))
    print("Dataset names:")
    print(names(datasets_geno))
    print(names(datasets_pheno))

    process_dataset <- function(genotype_df, phenotype_df, data_source) {
        print(paste("Processing dataset for:", data_source))
        crs_list_geno <- list()
        crs_list_pheno <- list()
        year_list <- list()

        # Get common years between genotype and phenotype datasets
        geno_years <- sort(unique(genotype_df$Year))
        pheno_years <- sort(unique(phenotype_df$Year))
        available_years <- sort(intersect(geno_years, pheno_years))
        
        print(paste("Available years:", paste(available_years, collapse = ", ")))

        for (year in available_years) {
            print(paste("Processing year:", year))
            i <- match(year, available_years)
            
            tryCatch({
                # Get rules
                if (rules_selected == "best") {
                    rules_pheno <- implement_cuts(df = phenotype_df, resistance_indicator = "phenotype", 
                                               target = target, cut_off = cut_off, 
                                               measures_used = measures_used, year = year, agg = FALSE)
                    rules_geno <- implement_cuts(df = genotype_df, resistance_indicator = "genotype", 
                                              target = target, cut_off = cut_off, 
                                              measures_used = measures_used, year = year, agg = FALSE)
                } else if (rules_selected == "all") {
                    rules_pheno <- get_rules_or_itemsets(df = phenotype_df, resistance_indicator = "phenotype", 
                                                       year = year, target = target, agg = FALSE)
                    rules_geno <- get_rules_or_itemsets(df = genotype_df, resistance_indicator = "genotype", 
                                                      year = year, target = target, agg = FALSE)
                }

                print(paste("Number of rules - Pheno:", length(rules_pheno), "Geno:", length(rules_geno)))

                # Only process if both rule sets have items
                if (length(rules_pheno) > 0 && length(rules_geno) > 0) {
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
                } else {
                    print(paste("Skipping year", year, "- insufficient rules"))
                }
            }, error = function(e) {
                print(paste("Error processing year", year, ":", e$message))
            })
        }

        # Only create dataframe if we have results
        if (length(year_list) > 0) {
            print("Creating data frame for dataset")
            result_df <- data.frame(
                year = c(unlist(year_list), unlist(year_list)),
                value = c(unlist(crs_list_geno), unlist(crs_list_pheno)),
                type = rep(c("Geno in Pheno", "Pheno in Geno"), each = length(year_list)),
                data_source = data_source
            )
            print(paste("Number of rows in result:", nrow(result_df)))
            return(result_df)
        } else {
            print("No valid data found for this dataset pair")
            return(NULL)
        }
    }

    all_crs_df <- data.frame()
    
    # Process each dataset pair
    for (dataset_name in names(datasets_geno)) {
        print(paste("Processing dataset pair:", dataset_name))
        base_name <- gsub(" Genotype", "", dataset_name)
        pheno_name <- paste(base_name, "Phenotype")
        
        crs_df <- process_dataset(
            genotype_df = datasets_geno[[dataset_name]],
            phenotype_df = datasets_pheno[[pheno_name]],
            data_source = base_name
        )
        
        all_crs_df <- rbind(all_crs_df, crs_df)
    }
    
    # Create the interaction variable for plotting
    all_crs_df$interaction <- paste(all_crs_df$data_source, all_crs_df$type)

 # Create visual mapping
   visual_mapping <- list(
       "Cecal Geno in Pheno" = list(color = "#0077BB", linetype = "solid", shape = 16, label = "Cecal Genotype in Phenotype"),
       "Cecal Pheno in Geno" = list(color = "#0077BB", linetype = "dashed", shape = 16, label = "Cecal Phenotype in Genotype"),
       "Retail Meats Geno in Pheno" = list(color = "#009E73", linetype = "solid", shape = 17, label = "Retail Meats Genotype in Phenotype"),
       "Retail Meats Pheno in Geno" = list(color = "#009E73", linetype = "dashed", shape = 17, label = "Retail Meats Phenotype in Genotype")
   )
    # Create vectors for scale_*_manual
   color_values <- sapply(visual_mapping, function(x) x$color)
   linetype_values <- sapply(visual_mapping, function(x) x$linetype)
   shape_values <- sapply(visual_mapping, function(x) x$shape)
   labels <- sapply(visual_mapping, function(x) x$label)
    plot1 <- ggplot(all_crs_df, aes(x = year, y = value, color = interaction, linetype = interaction, shape = interaction)) +
       geom_line(linewidth = 1.5) +
       geom_point(size = 4) +
       theme(
           text = element_text(size = 25),
           plot.title = element_text(hjust = 0.5),
           panel.background = element_rect(fill = 'white', colour = 'black'),
           legend.position = "bottom",
           legend.title = element_blank(),
           legend.text = element_text(size = 15),
           legend.key.size = unit(1.2, "cm"),
           legend.box = "vertical",
           legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
       ) +
       labs(y = "Cumulative Proportion Captured", x = "Year", 
            title = str_glue({title})) + 
       scale_x_continuous(breaks = round(seq(min(all_crs_df$year), max(all_crs_df$year), by = 2),1)) +
       scale_color_manual(
           values = color_values,
           labels = labels,
           breaks = names(visual_mapping)
       ) +
       scale_linetype_manual(
           values = linetype_values,
           labels = labels,
           breaks = names(visual_mapping)
       ) +
       scale_shape_manual(
           values = shape_values,
           labels = labels,
           breaks = names(visual_mapping)
       ) +
       ylim(0,1) +
       guides(
           color = guide_legend(nrow = 3, byrow = TRUE),
           linetype = guide_legend(nrow = 3, byrow = TRUE),
           shape = guide_legend(nrow = 3, byrow = TRUE)
       )

    ggsave(plot = plot1, 
           filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/between_indicators/c_proportion_captured/combined_{target}_captured_bidirectional_{rules_selected}.png"),
           width = 10, height = 8)
    
    return(plot1)
}




################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
#####################################################   BETWEEN DATA SOURCES      ##############################################################################################
################################################################################################################################################################################    
################################################################################################################################################################################


# Function to calculate and plot cumulative percent captured between cecal and retail for both genotype and phenotype
cecal_vs_retail_cumulative_percent_captured_multi <- function(cecal_genotype_df, cecal_phenotype_df, retail_genotype_df, retail_phenotype_df, class_level, target, rules_selected = "best", cut_off, measures_used, data_source1, data_source2, data_partition) {
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
    # Create combined plot with all four lines
    class_prefix <- if(class_level) "Class " else ""
    title_prefix <- if(class_level) "Class Level " else "Individual Antimicrobial/Gene "

    plot1 <- ggplot(crs_df) +
        geom_line(aes(x = year, y = d1in2, color = paste(indicator, data_source1, "in", data_source2)), linewidth = 2) +
        geom_point(aes(x = year, y = d1in2, color = paste(indicator, data_source1, "in", data_source2)), size = 3) +
        geom_line(aes(x = year, y = d2in1, color = paste(indicator, data_source2, "in", data_source1)), linewidth = 2) +
        geom_point(aes(x = year, y = d2in1, color = paste(indicator, data_source2, "in", data_source1)), size = 3) +
        theme(
            text = element_text(size = 25),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 15),
            legend.key.size = unit(0.8, "cm"),
            legend.box = "vertical",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
        ) +
        labs(y = "Cumulative Proportion Captured", x = "Year", 
             title = str_glue("{title_prefix}Cumulative Rule Set Comparison")) +
        scale_x_continuous(breaks = round(seq(min(crs_df$year), max(crs_df$year), by = 2),1)) +
        scale_color_manual(values = c("#1E90FF", "#009E73", "#F40808", "#FF8C00"),
                           labels = c(str_glue("Genotype Cecal in Retail Meats"), str_glue("Genotype Retail Meats in Cecal"), 
                                      str_glue("Phenotype Cecal in Retail Meats"), str_glue("Phenotype Retail Meats in Cecal"))) +
        ylim(0,1) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE))
    
    filename_suffix <- if(class_level == TRUE) "_CLASS" else ""
    ggsave(plot = plot1, filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/between_data_sources/c_proportion_captured/combined_{target}_captured_bidirectional_{rules_selected}{filename_suffix}.png"),
           width = 10, height = 8)
    return(plot1)
}







# Function to calculate and plot cumulative percent captured between cecal and retail for both genotype and phenotype
cecal_vs_retail_cumulative_percent_captured_multi_v2 <- function(datasets_geno, datasets_pheno, class_level, target, 
                                                               rules_selected = "best", cut_off, measures_used, data_partition) {
    process_dataset_pair <- function(df1, df2, indicator, data_source1, data_source2) {
        crs_list_1in2 <- list()
        crs_list_2in1 <- list()
        year_list <- list()

        # Use common years between the datasets
        common_years <- sort(intersect(unique(df1$Year), unique(df2$Year)))
        
        for (year in common_years) {
            i <- match(year, common_years)
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
            value = c(unlist(crs_list_1in2), unlist(crs_list_2in1)),
            type = rep(c(paste(data_source1, "in", data_source2), 
                        paste(data_source2, "in", data_source1)), 
                      each = length(year_list)),
            indicator = indicator
        )
    }

    all_crs_df <- data.frame()
    
    # Define all pairs of data sources to compare
    data_sources <- c("Cecal", "Retail Meats", "NAHLN")
    pairs <- combn(data_sources, 2, simplify = FALSE)
    
    # Process each pair for both genotype and phenotype
    for (pair in pairs) {
        source1 <- pair[1]
        source2 <- pair[2]
        
        # Get the corresponding datasets
        geno1_name <- paste(source1, "Genotype")
        geno2_name <- paste(source2, "Genotype")
        pheno1_name <- paste(source1, "Phenotype")
        pheno2_name <- paste(source2, "Phenotype")
        
        # Process genotype comparison
        geno_df <- process_dataset_pair(
            datasets_geno[[geno1_name]], 
            datasets_geno[[geno2_name]], 
            "genotype", source1, source2
        )
        
        # Process phenotype comparison
        pheno_df <- process_dataset_pair(
            datasets_pheno[[pheno1_name]], 
            datasets_pheno[[pheno2_name]], 
            "phenotype", source1, source2
        )
        
        all_crs_df <- rbind(all_crs_df, geno_df, pheno_df)
    }
    
    # Create the interaction variable for plotting
    all_crs_df$interaction <- paste(all_crs_df$indicator, all_crs_df$type)

    # Create color mapping
    color_mapping <- list(
        "genotype Cecal in Retail Meats" = list(color = "#1E90FF", label = "Genotype Cecal in Retail Meats"),
        "genotype Retail Meats in Cecal" = list(color = "#009E73", label = "Genotype Retail Meats in Cecal"),
        "genotype Cecal in NAHLN" = list(color = "#4169E1", label = "Genotype Cecal in NAHLN"),
        "genotype NAHLN in Cecal" = list(color = "#00CED1", label = "Genotype NAHLN in Cecal"),
        "genotype Retail Meats in NAHLN" = list(color = "#98FB98", label = "Genotype Retail Meats in NAHLN"),
        "genotype NAHLN in Retail Meats" = list(color = "#3CB371", label = "Genotype NAHLN in Retail Meats"),
        "phenotype Cecal in Retail Meats" = list(color = "#F40808", label = "Phenotype Cecal in Retail Meats"),
        "phenotype Retail Meats in Cecal" = list(color = "#FF8C00", label = "Phenotype Retail Meats in Cecal"),
        "phenotype Cecal in NAHLN" = list(color = "#FF4500", label = "Phenotype Cecal in NAHLN"),
        "phenotype NAHLN in Cecal" = list(color = "#FF69B4", label = "Phenotype NAHLN in Cecal"),
        "phenotype Retail Meats in NAHLN" = list(color = "#FFB6C1", label = "Phenotype Retail Meats in NAHLN"),
        "phenotype NAHLN in Retail Meats" = list(color = "#FF1493", label = "Phenotype NAHLN in Retail Meats")
    )

    # Create vectors for scale_color_manual
    color_values <- sapply(color_mapping, function(x) x$color)
    color_labels <- sapply(color_mapping, function(x) x$label)

    class_prefix <- if(class_level) "Class " else ""
    title_prefix <- if(class_level) "Class Level " else "Individual Antimicrobial/Gene "

    plot1 <- ggplot(all_crs_df, aes(x = year, y = value, color = interaction)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        theme(
            text = element_text(size = 25),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 15),
            legend.key.size = unit(0.8, "cm"),
            legend.box = "vertical",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
        ) +
        labs(y = "Cumulative Proportion Captured", x = "Year", 
             title = str_glue("{title_prefix}Cumulative Rule Set Comparison")) +
        scale_x_continuous(breaks = round(seq(min(all_crs_df$year), max(all_crs_df$year), by = 2),1)) +
        scale_color_manual(
            values = color_values,
            labels = color_labels,
            breaks = names(color_mapping)
        ) +
        ylim(0,1) +
        guides(color = guide_legend(nrow = 4, byrow = TRUE))

    filename_suffix <- if(class_level) "_CLASS" else ""
    ggsave(plot = plot1, 
           filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/between_data_sources/c_proportion_captured/combined_{target}_captured_bidirectional_{rules_selected}{filename_suffix}.png"),
           width = 10, height = 8)
    
    return(plot1)
}






# Function for genotype comparisons only
genotype_cumulative_percent_captured_multi <- function(datasets_geno, class_level, target, 
                                                     rules_selected = "best", cut_off, measures_used, data_partition) {
    process_dataset_pair <- function(df1, df2, data_source1, data_source2) {
        crs_list_1in2 <- list()
        crs_list_2in1 <- list()
        year_list <- list()

        # Use common years between the datasets
        common_years <- sort(intersect(unique(df1$Year), unique(df2$Year)))
        
        for (year in common_years) {
            i <- match(year, common_years)
            if (rules_selected == "best") {
                rules1 <- implement_cuts(df = df1, resistance_indicator = "genotype", target = target, 
                                      cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
                rules2 <- implement_cuts(df = df2, resistance_indicator = "genotype", target = target, 
                                      cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
            } else if (rules_selected == "all") {
                rules1 <- get_rules_or_itemsets(df = df1, resistance_indicator = "genotype", 
                                              year = year, target = target, agg = FALSE)
                rules2 <- get_rules_or_itemsets(df = df2, resistance_indicator = "genotype", 
                                              year = year, target = target, agg = FALSE)
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
            value = c(unlist(crs_list_1in2), unlist(crs_list_2in1)),
            type = rep(c(paste(data_source1, "in", data_source2), 
                        paste(data_source2, "in", data_source1)), 
                      each = length(year_list))
        )
    }

    all_crs_df <- data.frame()
    
    # Define all pairs of data sources to compare
    data_sources <- c("Cecal", "Retail Meats", "NAHLN")
    pairs <- combn(data_sources, 2, simplify = FALSE)
    
    # Process each pair for genotype only
    for (pair in pairs) {
        source1 <- pair[1]
        source2 <- pair[2]
        
        # Get the corresponding datasets
        geno1_name <- paste(source1, "Genotype")
        geno2_name <- paste(source2, "Genotype")
        
        # Process genotype comparison
        geno_df <- process_dataset_pair(
            datasets_geno[[geno1_name]], 
            datasets_geno[[geno2_name]], 
            source1, source2
        )
        
        all_crs_df <- rbind(all_crs_df, geno_df)
    }

    # Create visual mapping for the pairs
   visual_mapping <- list(
       "Cecal in Retail Meats" = list(color = "#0077BB", linetype = "solid", shape = 17),
       "Retail Meats in Cecal" = list(color = "#009E73", linetype = "solid", shape = 16),
       "Cecal in NAHLN" = list(color = "#0077BB", linetype = "solid", shape = 15),
       "NAHLN in Cecal" = list(color = "#EE7733", linetype = "dashed", shape = 16),
       "Retail Meats in NAHLN" = list(color = "#009E73", linetype = "dashed", shape = 15),
       "NAHLN in Retail Meats" = list(color = "#EE7733", linetype = "dashed", shape = 17)
   )
    # Create vectors for scale_*_manual
   color_values <- sapply(visual_mapping, function(x) x$color)
   linetype_values <- sapply(visual_mapping, function(x) x$linetype)
   shape_values <- sapply(visual_mapping, function(x) x$shape)
    class_prefix <- if(class_level) "Class " else ""
   title_prefix <- if(class_level) "Class Level " else "Individual Antimicrobial/Gene "
    
    
    # Add small jitter to y-values for better visibility of overlapping lines
    all_crs_df$value_jittered <- all_crs_df$value + runif(nrow(all_crs_df), -0.01, 0.01)

    plot1 <- ggplot(all_crs_df, aes(x = year, y = value_jittered, color = type, linetype = type, shape = type)) +
        geom_line(linewidth = 1.5, alpha = 0.7) +
        geom_point(size = 4) +
        theme(
            text = element_text(size = 25),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1.2, "cm"),
            legend.box = "vertical",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
        ) +
        labs(y = "Cumulative Proportion Captured", x = "Year", 
             title = str_glue("{title_prefix}Genotype Cumulative Rule Set Comparison")) +
        scale_x_continuous(breaks = round(seq(min(all_crs_df$year), max(all_crs_df$year), by = 2),1)) +
        scale_color_manual(values = color_values) +
        scale_linetype_manual(values = linetype_values) +
        scale_shape_manual(values = shape_values) +
        ylim(-0.05, 1.05) +  # Slightly expanded y-limits to accommodate jitter
        guides(
            color = guide_legend(nrow = 2, byrow = TRUE),
            linetype = guide_legend(nrow = 2, byrow = TRUE),
            shape = guide_legend(nrow = 2, byrow = TRUE)
        )

    filename_suffix <- if(class_level) "_CLASS" else ""
    ggsave(plot = plot1, 
           filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/between_data_sources/c_proportion_captured/combined_{target}_captured_bidirectional_{rules_selected}_genotype{filename_suffix}.png"),
           width = 10, height = 8)
    
    return(plot1)
}




# Function for phenotype comparisons only
phenotype_cumulative_percent_captured_multi <- function(datasets_pheno, class_level, target, 
                                                      rules_selected = "best", cut_off, measures_used, data_partition) {
    process_dataset_pair <- function(df1, df2, data_source1, data_source2) {
        crs_list_1in2 <- list()
        crs_list_2in1 <- list()
        year_list <- list()

        # Use common years between the datasets
        common_years <- sort(intersect(unique(df1$Year), unique(df2$Year)))
        
        for (year in common_years) {
            i <- match(year, common_years)
            if (rules_selected == "best") {
                rules1 <- implement_cuts(df = df1, resistance_indicator = "phenotype", target = target, 
                                      cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
                rules2 <- implement_cuts(df = df2, resistance_indicator = "phenotype", target = target, 
                                      cut_off = cut_off, measures_used = measures_used, year = year, agg = FALSE)
            } else if (rules_selected == "all") {
                rules1 <- get_rules_or_itemsets(df = df1, resistance_indicator = "phenotype", 
                                              year = year, target = target, agg = FALSE)
                rules2 <- get_rules_or_itemsets(df = df2, resistance_indicator = "phenotype", 
                                              year = year, target = target, agg = FALSE)
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
            value = c(unlist(crs_list_1in2), unlist(crs_list_2in1)),
            type = rep(c(paste(data_source1, "in", data_source2), 
                        paste(data_source2, "in", data_source1)), 
                      each = length(year_list))
        )
    }

    all_crs_df <- data.frame()
    
    # Define all pairs of data sources to compare
    data_sources <- c("Cecal", "Retail Meats", "NAHLN")
    pairs <- combn(data_sources, 2, simplify = FALSE)
    
    # Process each pair for phenotype only
    for (pair in pairs) {
        source1 <- pair[1]
        source2 <- pair[2]
        
        # Get the corresponding datasets
        pheno1_name <- paste(source1, "Phenotype")
        pheno2_name <- paste(source2, "Phenotype")
        
        # Process phenotype comparison
        pheno_df <- process_dataset_pair(
            datasets_pheno[[pheno1_name]], 
            datasets_pheno[[pheno2_name]], 
            source1, source2
        )
        
        all_crs_df <- rbind(all_crs_df, pheno_df)
    }

    # Create visual mapping for the pairs
   visual_mapping <- list(
       "Cecal in Retail Meats" = list(color = "#0077BB", linetype = "dashed", shape = 17),
       "Retail Meats in Cecal" = list(color = "#009E73", linetype = "dashed", shape = 16),
       "Cecal in NAHLN" = list(color = "#0077BB", linetype = "dashed", shape = 15),
       "NAHLN in Cecal" = list(color = "#EE7733", linetype = "solid", shape = 16),
       "Retail Meats in NAHLN" = list(color = "#009E73", linetype = "solid", shape = 15),
       "NAHLN in Retail Meats" = list(color = "#EE7733", linetype = "solid", shape = 17)
   )
    # Create vectors for scale_*_manual
   color_values <- sapply(visual_mapping, function(x) x$color)
   linetype_values <- sapply(visual_mapping, function(x) x$linetype)
   shape_values <- sapply(visual_mapping, function(x) x$shape)
    class_prefix <- if(class_level) "Class " else ""
   title_prefix <- if(class_level) "Class Level " else "Individual Antimicrobial/Gene "
    
    # Add small jitter to y-values for better visibility of overlapping lines
    all_crs_df$value_jittered <- all_crs_df$value + runif(nrow(all_crs_df), -0.01, 0.01)

    plot1 <- ggplot(all_crs_df, aes(x = year, y = value_jittered, color = type, linetype = type, shape = type)) +
        geom_line(linewidth = 1.5, alpha = 0.7) +
        geom_point(size = 4) +
        theme(
            text = element_text(size = 25),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1.2, "cm"),
            legend.box = "vertical",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
        ) +
        labs(y = "Cumulative Proportion Captured", x = "Year", 
             title = str_glue("{title_prefix}Phenotype Cumulative Rule Set Comparison")) +
        scale_x_continuous(breaks = round(seq(min(all_crs_df$year), max(all_crs_df$year), by = 2),1)) +
        scale_color_manual(values = color_values) +
        scale_linetype_manual(values = linetype_values) +
        scale_shape_manual(values = shape_values) +
        ylim(-0.05, 1.05) +  # Slightly expanded y-limits to accommodate jitter
        guides(
            color = guide_legend(nrow = 2, byrow = TRUE),
            linetype = guide_legend(nrow = 2, byrow = TRUE),
            shape = guide_legend(nrow = 2, byrow = TRUE)
        )


    filename_suffix <- if(class_level) "_CLASS" else ""
    ggsave(plot = plot1, 
           filename = str_glue("combined_percent_captured/{data_partition}/figures/shared_rules_analysis/between_data_sources/c_proportion_captured/combined_{target}_captured_bidirectional_{rules_selected}_phenotype{filename_suffix}.png"),
           width = 10, height = 8)
    
    return(plot1)
}





