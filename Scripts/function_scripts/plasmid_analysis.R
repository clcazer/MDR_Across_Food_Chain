library(stringr)
library(rlist)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(readxl)
library(cowplot)
library(arules)
library(dplyr)
library(reshape2)
library(gridExtra)
library(grid)
library(lubridate)
library(dplyr)
library(tidyr)
library(purrr)








comp_rules_to_plasmids <- function(rules_df, cecal_plasmids_df, retail_plasmids_df, nahln_plasmids_df) {
    # Add data source to each plasmid dataframe
    cecal_plasmids_df$Data_Source <- "Cecal"
    retail_plasmids_df$Data_Source <- "Retail_Meats"
    nahln_plasmids_df$Data_Source <- "NAHLN"
    
    # Function to normalize combinations
    normalize_combination <- function(combination) {
        combination <- gsub("[{}'\"]", "", combination)
        genes <- unlist(strsplit(combination, "[,|]"))
        genes <- sort(trimws(genes))
        genes <- genes[genes != ""]
        return(paste(genes, collapse = ","))
    }
    
    # Function to deduplicate plasmid dataframe
    deduplicate_plasmids <- function(df) {
        # Add normalized combinations
        df$normalized_combination <- sapply(df$combination, normalize_combination)
        
        # Group by normalized combination and sum contigs
        deduped <- aggregate(
            num_contigs_with_combo ~ normalized_combination + Data_Source, 
            data = df,
            FUN = sum
        )
        
        # Keep the original combination format (take first occurrence)
        combinations <- df[!duplicated(df$normalized_combination), c("combination", "normalized_combination")]
        deduped <- merge(deduped, combinations, by = "normalized_combination")
        
        return(deduped)
    }
    
    # Deduplicate each plasmid dataset
    cecal_plasmids_df <- deduplicate_plasmids(cecal_plasmids_df)
    retail_plasmids_df <- deduplicate_plasmids(retail_plasmids_df)
    nahln_plasmids_df <- deduplicate_plasmids(nahln_plasmids_df)
    
    # Combine the deduplicated dataframes
    plasmids_df <- rbind(cecal_plasmids_df, retail_plasmids_df, nahln_plasmids_df)
    
    # Normalize rules dataframe
    rules_df$normalized_combination <- sapply(rules_df$gene_combination, normalize_combination)
    
    # Create a result data frame with pre-allocated columns
    results <- data.frame(
        Data_Source = character(),
        Year = numeric(),
        rule_combination = character(),
        class_combination = character(),
        plasmid_combination = character(),
        support = numeric(),
        confidence = numeric(),
        lift = numeric(),
        num_contigs = numeric(),
        mismatch_distance = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Compare combinations, matching only within the same data source
    for (i in 1:nrow(rules_df)) {
        rule_combination <- rules_df$normalized_combination[i]
        rule_genes <- unlist(strsplit(rule_combination, ","))
        
        # Filter plasmids_df to only include rows from the same data source
        matching_plasmids <- plasmids_df[plasmids_df$Data_Source == rules_df$Data_Source[i], ]
        
        for (j in 1:nrow(matching_plasmids)) {
            plasmid_combination <- matching_plasmids$normalized_combination[j]
            plasmid_genes <- unlist(strsplit(plasmid_combination, ","))
            
            # Calculate Jaccard similarity
            intersection <- length(intersect(rule_genes, plasmid_genes))
            union <- length(unique(c(rule_genes, plasmid_genes)))
            mismatch_distance <- 1 - (intersection / union)
            
            if (mismatch_distance <= 0.7) {
                new_row <- data.frame(
                    Data_Source = rules_df$Data_Source[i],
                    Year = rules_df$Year[i],
                    rule_combination = rules_df$gene_combination[i],
                    class_combination = rules_df$class_combination[i],
                    plasmid_combination = matching_plasmids$combination[j],
                    support = rules_df$support[i],
                    confidence = rules_df$confidence[i],
                    lift = rules_df$lift[i],
                    num_contigs = matching_plasmids$num_contigs_with_combo[j],
                    mismatch_distance = mismatch_distance,
                    stringsAsFactors = FALSE
                )
                results <- rbind(results, new_row)
            }
        }
    }
    
    # Keep only the best match for each plasmid combination within each data source
    if (nrow(results) > 0) {
        # Create two copies of the results for different sorting
        plasmid_best_matches <- results
        rule_best_matches <- results
        
        # Best matches for each plasmid combination
        plasmid_best_matches <- plasmid_best_matches[order(plasmid_best_matches$Data_Source, 
                               plasmid_best_matches$plasmid_combination,
                               plasmid_best_matches$mismatch_distance,
                               -plasmid_best_matches$support,
                               -plasmid_best_matches$confidence,
                               -plasmid_best_matches$lift), ]
        
        plasmid_best_matches <- plasmid_best_matches[!duplicated(plasmid_best_matches[c("Data_Source", "plasmid_combination")]), ]
        
        # Best matches for each rule combination
        rule_best_matches <- rule_best_matches[order(rule_best_matches$Data_Source, 
                               rule_best_matches$rule_combination,
                               rule_best_matches$mismatch_distance,
                               -rule_best_matches$num_contigs), ]
        
        rule_best_matches <- rule_best_matches[!duplicated(rule_best_matches[c("Data_Source", "rule_combination")]), ]
        
        # Write both sets of results to separate files
        write.csv(plasmid_best_matches, "Virulence_Plasmid_outputs/plasmid_rule_comp/plasmid_rule_comp_plasmid_best.csv", row.names = FALSE)
        write.csv(rule_best_matches, "Virulence_Plasmid_outputs/plasmid_rule_comp/plasmid_rule_comp_rule_best.csv", row.names = FALSE)
        
        # Return both sets of results as a list
        return(list(
            plasmid_best_matches = plasmid_best_matches,
            rule_best_matches = rule_best_matches
        ))
    } else {
        return(data.frame(message = "No matches found"))
    }
}


rules_df <- read.csv("combined_csv_outputs/Unique_top_1000_rules_by_dataset_with_classes.csv")
cecal_plasmids <- read.csv("Raw_Data/Virulence_Plasmid_data/plasmid-amr-NARMS_CECAL.csv")
retail_plasmids <- read.csv("Raw_Data/Virulence_Plasmid_data/plasmid-amr-NARMS_MEAT.csv")
nahln_plasmids <- read.csv("Raw_Data/Virulence_Plasmid_data/plasmid-amr-NAHLN.csv")

plasmid_comp_res <- comp_rules_to_plasmids(rules_df, cecal_plasmids, retail_plasmids, nahln_plasmids)



plasmid_rule_match_summary <- function(matched_df, match_type = c("plasmid", "rule")) {
    # Validate and process match_type parameter
    match_type <- match.arg(match_type)
    file_suffix <- ifelse(match_type == "plasmid", "plasmid_best", "rule_best")
    
    # Basic summaries by Data_Source
    basic_summaries <- matched_df %>%
        group_by(Data_Source) %>%
        summarize(
            n_matches = n(),
            n_unique_plasmids = n_distinct(plasmid_combination),
            n_unique_rules = n_distinct(rule_combination),
            
            # Mismatch distance statistics
            avg_mismatch = mean(mismatch_distance, na.rm = TRUE),
            sd_mismatch = sd(mismatch_distance, na.rm = TRUE),
            
            # Correlation between num_contigs and support
            contig_support_cor = cor(num_contigs, support, method = "spearman", 
                                   use = "complete.obs"),
            
            # Basic metrics
            avg_support = mean(support, na.rm = TRUE),
            sd_support = sd(support, na.rm = TRUE),
            avg_confidence = mean(confidence, na.rm = TRUE),
            avg_lift = mean(lift, na.rm = TRUE),
            avg_contigs = mean(num_contigs, na.rm = TRUE),
            pct_perfect_match = mean(mismatch_distance == 0, na.rm = TRUE) * 100
        ) %>%
        mutate(across(where(is.numeric), ~round(., 3)))

    # Function to convert string combinations to character vectors
    parse_combination <- function(combo_str) {
        gsub("[{}]", "", combo_str) %>%
            strsplit(",") %>%
            unlist() %>%
            trimws()
    }

    # Analyze mismatch types
    mismatch_types <- matched_df %>%
        mutate(
            rule_genes = map(rule_combination, parse_combination),
            plasmid_genes = map(plasmid_combination, ~strsplit(., "\\|")[[1]]),
            mismatch_type = case_when(
                mismatch_distance == 0 ~ "perfect_match",
                map2_lgl(plasmid_genes, rule_genes, 
                        ~all(.x %in% .y)) ~ "plasmid_subset_of_rule",
                map2_lgl(rule_genes, plasmid_genes, 
                        ~all(.x %in% .y)) ~ "rule_subset_of_plasmid",
                TRUE ~ "pure_mismatch"
            )
        ) %>%
        group_by(Data_Source) %>%
        count(mismatch_type) %>%
        pivot_wider(names_from = mismatch_type, 
                   values_from = n, 
                   values_fill = 0)

    # Class-level analysis
    class_summaries <- matched_df %>%
        group_by(Data_Source, class_combination) %>%
        summarize(
            n = n(),
            avg_support = mean(support),
            avg_confidence = mean(confidence),
            avg_mismatch = mean(mismatch_distance),
            .groups = "drop"
        ) %>%
        arrange(Data_Source, desc(n)) %>%
        group_by(Data_Source) %>%
        slice_head(n = 5)  # Top 5 class combinations per source

    # Gene-level analysis
    gene_summaries <- matched_df %>%
        mutate(rule_genes = map(rule_combination, parse_combination)) %>%
        unnest(rule_genes) %>%
        group_by(Data_Source, rule_genes) %>%
        summarize(
            frequency = n(),
            avg_support = mean(support),
            avg_confidence = mean(confidence),
            .groups = "drop"
        ) %>%
        arrange(Data_Source, desc(frequency)) %>%
        group_by(Data_Source) %>%
        slice_head(n = 10)  # Top 10 genes per source

    # Save all summaries to separate CSV files with appropriate suffixes
    write.csv(basic_summaries, 
              sprintf("Virulence_Plasmid_outputs/plasmid_rule_comp/%s_basic_summary.csv", file_suffix), 
              row.names = FALSE)
    write.csv(mismatch_types, 
              sprintf("Virulence_Plasmid_outputs/plasmid_rule_comp/%s_mismatch_types.csv", file_suffix), 
              row.names = FALSE)
    write.csv(class_summaries, 
              sprintf("Virulence_Plasmid_outputs/plasmid_rule_comp/%s_class_summary.csv", file_suffix), 
              row.names = FALSE)
    write.csv(gene_summaries, 
              sprintf("Virulence_Plasmid_outputs/plasmid_rule_comp/%s_gene_summary.csv", file_suffix), 
              row.names = FALSE)

    # Return list of all summaries
    return(list(
        basic_summaries = basic_summaries,
        mismatch_types = mismatch_types,
        class_summaries = class_summaries,
        gene_summaries = gene_summaries
    ))
}


best_plasmid_matches <- read.csv("Virulence_Plasmid_outputs/plasmid_rule_comp/plasmid_rule_comp_plasmid_best.csv")
best_rule_matches <- read.csv("Virulence_Plasmid_outputs/plasmid_rule_comp/plasmid_rule_comp_rule_best.csv")

plasmid_rule_match_summary(matched_df = best_plasmid_matches, match_type = "plasmid")
plasmid_rule_match_summary(matched_df = best_rule_matches, match_type = "rule")