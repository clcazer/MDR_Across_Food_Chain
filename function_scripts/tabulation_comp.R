library(arules)
library(here)
library(stringr)


source(here("function_scripts", "data_wrangling.R"))
source(here("function_scripts", "rule_mining_and_selection.R"))
source(here("function_scripts", "rule_analysis.R"))


tabulate_and_compare <- function(df, rules_selected, cut_off, measures_used, resistance_indicator, data_source, class_level) {
     # Remove the MER column before processing
   # df <- df[, names(df) != "MER"]
    #df <- df[, names(df) != "SMX"]
    #df <- df[, names(df) != "AZI"]

    #gsub BETA.LACTAM with BETA-LACTAM in df colnames
    colnames(df) <- gsub("BETA.LACTAM", "BETA-LACTAM", colnames(df))
    #gsub FOLATE.PATHWAY.INHIBITOR with FOLATE-PATHWAY-INHIBITOR in df colnames
    colnames(df) <- gsub("FOLATE.PATHWAY.INHIBITOR", "FOLATE-PATHWAY-INHIBITOR", colnames(df))


    gene_class_mappings <- read.csv("gene_class_mappings.csv")


    
    # Create a copy of df for tabulation purposes only
    df_tabulation <- df
    df_tabulation[is.na(df_tabulation)] <- 0
    
    # Create a new column for resistance patterns using df_tabulation
    df_tabulation$resistance_pattern <- apply(df_tabulation, 1, function(row) {
        # Get the names of the drugs the isolate is resistant to
        resistant_drugs <- names(df_tabulation)[which(row == 1)]
        # Sort drug names for consistency and combine them into a single string
        paste(sort(resistant_drugs), collapse = ", ")
    })
    
    # Count the frequency of each unique resistance pattern
    resistance_summary <- as.data.frame(table(df_tabulation$resistance_pattern), stringsAsFactors = FALSE)
    colnames(resistance_summary) <- c("Resistance_Pattern", "Frequency")
    
    # Order by frequency (descending)
    resistance_summary <- resistance_summary[order(-resistance_summary$Frequency), ]
    
    # Reset row names
    rownames(resistance_summary) <- NULL


    # Initialize empty list to store rule dataframes from each year
    all_rules_df <- list()
    
    # Add data validation
    # Get only numeric columns, excluding ID and Year
    numeric_cols <- names(df_tabulation)[!names(df_tabulation) %in% c("ID", "Year", "resistance_pattern")]
    numeric_data <- as.matrix(df_tabulation[, numeric_cols])

    # Add more detailed debugging
    print("Dimensions of numeric data:")
    print(dim(numeric_data))
    
    print("Class of numeric data:")
    print(class(numeric_data))
    
    print("Summary of numeric data:")
    print(summary(as.vector(numeric_data)))
    
    # Check for specific problematic values
    if (any(is.na(numeric_data))) {
        stop("Data contains NA values")
    }
    if (any(is.infinite(numeric_data))) {
        stop("Data contains infinite values")
    }
    if (any(!numeric_data %in% c(0,1))) {
        stop("Data contains values other than 0 and 1")
    }


    
    # Process rules year by year
    # Replace min/max with unique years to handle gaps
    unique_years <- sort(unique(df$Year))
    for (year in unique_years) {
        tryCatch({
            if (rules_selected == "best") {
                rules <- implement_cuts(df = df, resistance_indicator = resistance_indicator, 
                                      target = "rules", cut_off = cut_off, year = year, 
                                      agg = FALSE, measures_used = measures_used)
            } else if (rules_selected == "all") {
                rules <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, 
                                             year = year, target = "rules", agg = FALSE)
            }
            
            # If rules exist for this year, process them into a dataframe
            if (!is.null(rules) && class(rules) == "rules") {
                # Get LHS and RHS
                lhs_items <- labels(lhs(rules))
                rhs_items <- labels(rhs(rules))
                
                # Get quality measures for all rules at once
                quality_measures <- quality(rules)
                
                year_df <- data.frame(
                    lhs = lhs_items,
                    rhs = rhs_items,
                    cosine = quality_measures$cosine,
                    jaccard = quality_measures$jaccard,
                    kulczynski = quality_measures$kulczynski,
                    support = quality_measures$support,
                    stringsAsFactors = FALSE
                )
                
                # Add complete patterns for this year
                for (i in 1:nrow(year_df)) {
                    # Remove "{" and "}" from the strings and split
                    lhs_set <- strsplit(gsub("[{}]", "", year_df$lhs[i]), ",")[[1]]
                    rhs_set <- strsplit(gsub("[{}]", "", year_df$rhs[i]), ",")[[1]]
                    # Combine items, remove whitespace, and sort
                    all_items <- unique(trimws(c(lhs_set, rhs_set)))
                    year_df$complete_pattern[i] <- paste(sort(all_items), collapse = ", ")
                }
                
                all_rules_df[[as.character(year)]] <- year_df
            }
        }, error = function(e) {
            warning(sprintf("Error processing year %d: %s", year, e$message))
            return(NULL)
        })
    }
    
    # Combine all rule dataframes
    if (length(all_rules_df) == 0) {
        warning("No valid rules found")
        resistance_summary$recovered_resistance_pattern <- "No"
        return(resistance_summary)
    }
    
    # Combine all years' dataframes
    combined_rules_df <- do.call(rbind, all_rules_df)
    
    # Add recovered_resistance_pattern column to resistance_summary
    resistance_summary$recovered_resistance_pattern <- sapply(
        as.character(resistance_summary$Resistance_Pattern),
        function(pattern) {
            # Convert pattern to set of drugs
            pattern_drugs <- strsplit(pattern, ", ")[[1]]
            pattern_drug_set <- sort(pattern_drugs)
            
            # If only one drug, return NA
            if (length(pattern_drug_set) == 1) {
                return(NA)
            }
            
            # Check each rule pattern for subset relationship
            for (rule_pattern in combined_rules_df$complete_pattern) {
                rule_drugs <- strsplit(rule_pattern, ", ")[[1]]
                rule_drug_set <- sort(rule_drugs)
                
                # Check if either set is a subset of the other
                if (all(rule_drug_set %in% pattern_drug_set) #|| 
                    # all(pattern_drug_set %in% rule_drug_set)
                    ) {
                    return("Yes")
                }
            }
            return("No")
        }
    )
    
    # Before writing to CSV, clean up gene names if needed
    if (resistance_indicator == "genotype" && class_level == FALSE) {
        # Create a lookup table for gene name conversion
        gene_lookup <- setNames(
            gene_class_mappings$Gene_family,
            gsub("[^[:alnum:]]", ".", gene_class_mappings$Gene_family)
        )
        
        # Function to replace gene names in a pattern string
        replace_gene_names <- function(pattern) {
            if (is.na(pattern)) return(pattern)
            
            # Split pattern into individual genes
            genes <- strsplit(pattern, ", ")[[1]]
            
            # Replace each gene name with its original format
            original_genes <- sapply(genes, function(gene) {
                # Look up the original name
                if (gene %in% names(gene_lookup)) {
                    return(gene_lookup[gene])
                }
                return(gene)  # If no match found, return original
            })
            
            # Ensure we're working with character vector before sorting
            original_genes <- unlist(original_genes)
            
            # Recombine into pattern
            paste(sort(original_genes), collapse = ", ")
        }
        
        # Apply the replacement to both columns
        resistance_summary$Resistance_Pattern <- sapply(resistance_summary$Resistance_Pattern, replace_gene_names)
        combined_rules_df$complete_pattern <- sapply(combined_rules_df$complete_pattern, replace_gene_names)
    }
    
    # Add columns for quality measures
    resistance_summary$avg_cosine <- NA
    resistance_summary$avg_jaccard <- NA
    resistance_summary$avg_kulczynski <- NA
    resistance_summary$avg_support <- NA
    
    # Calculate average quality measures
    for (i in 1:nrow(resistance_summary)) {
        if (!is.na(resistance_summary$recovered_resistance_pattern[i]) && 
            resistance_summary$recovered_resistance_pattern[i] == "Yes") {
            pattern <- resistance_summary$Resistance_Pattern[i]
            pattern_drugs <- strsplit(pattern, ", ")[[1]]
            
            # Find all matching rules across all years
            matching_measures <- list(cosine = numeric(), 
                                    jaccard = numeric(),
                                    kulczynski = numeric(),
                                    support = numeric())
            
            for (year_df in all_rules_df) {
                if (!is.null(year_df)) {
                    for (j in 1:nrow(year_df)) {
                        rule_drugs <- strsplit(year_df$complete_pattern[j], ", ")[[1]]
                        if (all(rule_drugs %in% pattern_drugs)) {
                            matching_measures$cosine <- c(matching_measures$cosine, year_df$cosine[j])
                            matching_measures$jaccard <- c(matching_measures$jaccard, year_df$jaccard[j])
                            matching_measures$kulczynski <- c(matching_measures$kulczynski, year_df$kulczynski[j])
                            matching_measures$support <- c(matching_measures$support, year_df$support[j])
                        }
                    }
                }
            }
            
            # Calculate averages
            resistance_summary$avg_cosine[i] <- mean(matching_measures$cosine, na.rm = TRUE)
            resistance_summary$avg_jaccard[i] <- mean(matching_measures$jaccard, na.rm = TRUE)
            resistance_summary$avg_kulczynski[i] <- mean(matching_measures$kulczynski, na.rm = TRUE)
            resistance_summary$avg_support[i] <- mean(matching_measures$support, na.rm = TRUE)
        }
    }

    print(head(resistance_summary, 20))

    filename_suffix <- ifelse(class_level == TRUE, "_CLASS", "")
    write.csv(resistance_summary, str_glue("{data_source}/tabulation_comp/{resistance_indicator}/{data_source}_{resistance_indicator}_tabulation{filename_suffix}.csv"), row.names = FALSE)

    return(resistance_summary)

}



# get all the data
#retail
Retail_Meats_phenotype_df <- read.csv('Retail_Meats/NEW_Retail_Meats_wide_resStatus_phenotype.csv')

#cecal
cecal_phenotype_df <- read.csv("cecal/NEW_cecal_wide_resStatus_phenotype.csv")

#NAHLN
NAHLN_phenotype_df <- read.csv("NAHLN/NAHLN_wide_resStatus_phenotype.csv")




#run tabulation function for each df
retail_tab_comp <- tabulate_and_compare(df = Retail_Meats_phenotype_df, 
rules_selected = "best", cut_off = c(0.5, 0, 0.5, 0), 
measures_used = c("cosine", "jaccard", "kulczynski", "support"),
resistance_indicator = "phenotype",
data_source = "Retail_Meats",
class_level = FALSE)









cecal_tab_comp <- tabulate_and_compare(df = cecal_phenotype_df, 
rules_selected = "best", cut_off = c(0.5, 0, 0.5, 0), 
measures_used = c("cosine", "jaccard", "kulczynski", "support"),
resistance_indicator = "phenotype",
data_source = "cecal",
class_level = FALSE)









NAHLN_tab_comp <- tabulate_and_compare(df = NAHLN_phenotype_df, 
rules_selected = "best", cut_off = c(0.5, 0, 0.5, 0), 
measures_used = c("cosine", "jaccard", "kulczynski", "support"),
resistance_indicator = "genotype",
data_source = "NAHLN",
class_level = FALSE)







# Combine all three dataframes
# Add a 'Source' column to identify the origin of each row
retail_tab_comp$Source <- "Retail Meats"
cecal_tab_comp$Source <- "Cecal"
NAHLN_tab_comp$Source <- "NAHLN"

# Combine the dataframes
combined_df <- rbind(retail_tab_comp, cecal_tab_comp, NAHLN_tab_comp)

# Filter out patterns with less than two drugs
combined_df <- combined_df[sapply(strsplit(combined_df$Resistance_Pattern, ", "), length) >= 2, ]


# Sort by Frequency (descending) and keep top 20 rows
combined_df <- combined_df[order(-combined_df$Frequency), ][1:20, ]

# Round the average columns to 2 decimal places
combined_df$avg_cosine <- round(combined_df$avg_cosine, 2)
combined_df$avg_jaccard <- round(combined_df$avg_jaccard, 2)
combined_df$avg_kulczynski <- round(combined_df$avg_kulczynski, 2)
combined_df$avg_support <- round(combined_df$avg_support, 2)


write.csv(combined_df, "combined_tabComp_top20.csv", row.names = FALSE)
