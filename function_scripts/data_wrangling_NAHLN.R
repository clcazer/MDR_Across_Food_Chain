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





#######################################################
##########NAHLN PHENOTYPE DATA WRANGLING###############
#######################################################

wrangle_NAHLN_phenotype <- function(df, break_points_df) {

    #replace spaces in colnames with periods
    colnames(df) <- gsub(" ", ".", colnames(df))
    
    # Rename columns based on break_points_df mapping
    for(i in 1:nrow(break_points_df)) {
        full_name <- break_points_df$fullName[i]
        abbrev <- break_points_df$abbreviation[i]
        # Replace the full name with abbreviation in column names
        colnames(df) <- gsub(full_name, abbrev, colnames(df), fixed = TRUE)
    }

    # Clean MIC columns
    mic_cols <- grep("\\.MIC$", colnames(df))
    for(col in mic_cols) {
        # Convert to character first to handle any non-numeric entries
        df[[col]] <- as.character(df[[col]])
        # Extract only the numeric values (including decimals)
        df[[col]] <- as.numeric(gsub("[^0-9.]", "", df[[col]]))
    }
    
    # Clean Interpretation columns
    interp_cols <- grep("\\.Interpretation$", colnames(df))
    for(col in interp_cols) {
        df[[col]] <- case_when(
            df[[col]] == "R" ~ 1,
            df[[col]] == "S" ~ 0,
            df[[col]] == "NoCLSI" ~ NA_real_,
            TRUE ~ NA_real_  # Handle any other unexpected values
        )
    }



    # Select only desired columns (UniqueID, TestDate, and the MIC/Interpretation columns)
    desired_cols <- c("UniqueID", "TestDate", 
                     colnames(df)[grep("\\.(MIC|Interpretation)$", colnames(df))])
    df <- df[, desired_cols]

    # Convert TestDate to Year
    df$Year <- lubridate::year(as.Date(df$TestDate))
    
    # Remove TestDate column
    df <- df[, !colnames(df) %in% "TestDate"]
    #gsub TYL.tartrate column name to just TYL
    colnames(df) <- gsub("TYL\\.tartrate", "TYL", colnames(df))


    # Create new resistance incidence dataframe
    # Extract drug names from MIC columns (removing the .MIC suffix)
    drug_names <- gsub("\\.MIC$", "", colnames(df)[grep("\\.MIC$", colnames(df))])
    
    # Create empty resistance dataframe
    resistance_df <- data.frame(
        ID = df$UniqueID,
        Year = df$Year
    )
    
    # Add drug columns initialized with NA
    for(drug in drug_names) {
        resistance_df[[drug]] <- NA_real_
    }
    print(colnames(resistance_df))

    
    # Fill in resistance values
    for(drug in drug_names) {
        # Find if drug has a breakpoint
        breakpoint_row <- which(break_points_df$abbreviation == drug)
        
        if(length(breakpoint_row) > 0 && !is.na(break_points_df$resistant_value[breakpoint_row])) {
            # Use MIC comparison if breakpoint exists
            resistant_value <- break_points_df$resistant_value[breakpoint_row]
            mic_col <- paste0(drug, ".MIC")
            
            resistance_df[[drug]] <- ifelse(
                !is.na(df[[mic_col]]),
                ifelse(df[[mic_col]] >= resistant_value, 1, 0),
                NA
            )
        } else {
            # Use interpretation if no breakpoint
            interp_col <- paste0(drug, ".Interpretation")
            resistance_df[[drug]] <- df[[interp_col]]
        }
    }
    write.csv(resistance_df, "NAHLN/NAHLN_wide_resStatus_phenotype.csv", row.names =FALSE)
    return(resistance_df)
}







#######################################################
##########NAHLN GENOTYPE DATA WRANGLING################
#######################################################





keep_top_genes <- function(df, top_genes_df) {

    # Get top 50 genes from top_genes_df
    top_50_genes <- top_genes_df$Gene[1:30]
    
    # Get columns to keep (ID, Year, and top 50 genes)
    cols_to_keep <- c("UniqueID", "ID", "Year", top_50_genes)
    
    # Keep only those columns that exist in df
    existing_cols <- intersect(cols_to_keep, colnames(df))
    
    # Subset dataframe
    df <- df[, existing_cols]
    write.csv(df, "NAHLN/NAHLN_top30_family_level_genotype.csv")
    return(df)

}




NAHLN_process_genotype_data <- function(gene_df, Year_lookup_df, ID_lookup_df, data_source, output_level = 'gene') {
    gene_df <- gene_df[gene_df$meta == "NAHLN", ]
        #replace spaces in colnames with periods
    colnames(ID_lookup_df) <- gsub(" ", ".", colnames(ID_lookup_df))
    

    print("number of unique SAMNS")
    print(length(unique(gene_df$SAMN)))
    # Debug print for Year_lookup_df
    print("Sample of Year_lookup_df:")
    print(head(Year_lookup_df[, c("UniqueID", "TestDate")]))
    

    colnames(ID_lookup_df) <- gsub(" ", ".", colnames(ID_lookup_df))
    gene_df <- merge(gene_df,
                    ID_lookup_df[, c("BioSample", "Isolate.identifiers")],
                    by.x = "SAMN",
                    by.y = "BioSample",
                    all.x = TRUE)
    


    
   # Modified find_matching_id function to handle EC-Cow- prefix
   find_matching_id <- function(isolate_ids, unique_ids) {
       # Split isolate_ids into individual identifiers
       ids <- strsplit(isolate_ids, ",")[[1]]
       # Look specifically for pattern containing PPY
       ppy_ids <- ids[grep("PPY", ids)]
       # Clean up any quotes, whitespace, and remove EC-Cow- prefix
       ppy_ids <- gsub('"', '', trimws(ppy_ids))
       ppy_ids <- gsub("EC-Cow-", "", ppy_ids)
       # Return matching ID if found
       matching_id <- unique_ids[unique_ids %in% ppy_ids]
       if(length(matching_id) > 0) return(matching_id[1])
       return(NA)
   }
    
    
    # Convert TestDate to Year with error checking
    Year_lookup_df$Year <- tryCatch({
        year(as.Date(Year_lookup_df$TestDate, format = "%m/%d/%Y"))
    }, error = function(e) {
        print("Error converting dates. Trying alternative format...")
        year(as.Date(Year_lookup_df$TestDate, format = "%Y-%m-%d"))
    })
    
    # Create matched_unique_id column
    gene_df$matched_unique_id <- sapply(gene_df$Isolate.identifiers, 
                                      find_matching_id, 
                                      unique_ids = Year_lookup_df$UniqueID)

    
    # Try cleaning both columns before merge
    gene_df$matched_unique_id <- trimws(gene_df$matched_unique_id)
    Year_lookup_df$UniqueID <- trimws(Year_lookup_df$UniqueID)
    
    # Standardize encoding for both columns
    gene_df$matched_unique_id <- iconv(gene_df$matched_unique_id, to = "ASCII//TRANSLIT")
    Year_lookup_df$UniqueID <- iconv(Year_lookup_df$UniqueID, to = "ASCII//TRANSLIT")
    
    # Debug: Check a specific ID that should match
    print("Sample comparison:")
    sample_id <- gene_df$matched_unique_id[1]
    print(paste("Sample ID from gene_df:", sample_id))
    print("Matching rows in Year_lookup_df:")
    print(Year_lookup_df[Year_lookup_df$UniqueID == sample_id, ])
    
    # Create Year column from TestDate before merge
    Year_lookup_df$Year <- as.numeric(format(as.Date(Year_lookup_df$TestDate), "%Y"))
    

    
    # Merge with Year_lookup_df
    gene_df <- merge(gene_df,
                    Year_lookup_df[, c("UniqueID", "Year")],
                    by.x = "matched_unique_id",
                    by.y = "UniqueID",
                    all.x = TRUE)
    

    

    
    # Use the matched UniqueID as the final ID
    gene_df$ID <- gene_df$matched_unique_id
    
    # Print diagnostic information
    print(paste("Number of successful ID matches:", sum(!is.na(gene_df$matched_unique_id))))
    print(paste("Number of rows with ID/Year:", sum(!is.na(gene_df$ID) & !is.na(gene_df$Year))))



if (output_level == "gene") {
    # Remove duplicate gene entries for each isolate
    gene_df <- gene_df %>%
        distinct(ID, Element.symbol, .keep_all = TRUE)
    
    # Create the incidence matrix
    # First create a presence indicator
    gene_df$present <- 1
    

    
    # Pivot wider to create the incidence matrix
    incidence_matrix <- pivot_wider(gene_df,
                                  id_cols = c(ID, Year),
                                  names_from = Element.symbol,
                                  values_from = present,
                                  values_fill = 0)

    } else if (output_level == "class") {
    # Split entries in the Class column that contain multiple classes separated by /
    gene_df <- gene_df %>%
        separate_rows(Class, sep = "/")


        # Remove duplicate gene entries for each isolate
    gene_df <- gene_df %>%
        distinct(ID, Class, .keep_all = TRUE)


    # Create the incidence matrix
    # First create a presence indicator
    gene_df$present <- 1
    

    
    # Pivot wider to create the incidence matrix
    incidence_matrix <- pivot_wider(gene_df,
                                  id_cols = c(ID, Year),
                                  names_from = Class,
                                  values_from = present,
                                  values_fill = 0)

    } else if (output_level == "subclass") {
    # Split entries in the Class column that contain multiple classes separated by /
    gene_df <- gene_df %>%
        separate_rows(Subclass, sep = "/")


        # Remove duplicate gene entries for each isolate
    gene_df <- gene_df %>%
        distinct(ID, Subclass, .keep_all = TRUE)


    # Create the incidence matrix
    # First create a presence indicator
    gene_df$present <- 1
    

    
    # Pivot wider to create the incidence matrix
    incidence_matrix <- pivot_wider(gene_df,
                                  id_cols = c(ID, Year),
                                  names_from = Subclass,
                                  values_from = present,
                                  values_fill = 0)

    }

    
    incidence_matrix <- incidence_matrix %>%
        drop_na(ID, Year)


    # # Drop any column that has "cya" in it
    # incidence_matrix <- incidence_matrix %>%
    #     select(-matches("cya"))


    # Reorder columns to ensure ID and Year are first, followed by the rest in alphabetical order
    incidence_matrix <- incidence_matrix %>%
        select(ID, Year, sort(names(incidence_matrix)[-(1:2)]))

    if (output_level == "gene") {
    
    write.csv(incidence_matrix, file = str_glue("{data_source}/{data_source}_wide_corrected_genotype.csv"), row.names = FALSE)
    
    } else if (output_level == "class") {
    
    write.csv(incidence_matrix, file = str_glue("{data_source}/{data_source}_wide_class_level_genotype.csv"), row.names = FALSE)
    
    } else if (output_level == "subclass") {
    
    write.csv(incidence_matrix, file = str_glue("{data_source}/{data_source}_wide_Subclass_level_genotype.csv"), row.names = FALSE)
    
    }
    
    return(incidence_matrix)
}


# gene_df <- read.csv("Gene_Data/amrfinder_results-AMR-CORE.point-muts-included.csv")
# ID_lookup_df <- read_excel("NAHLN_Data/Isolate+BioSample.xlsx")
# Year_lookup_df <- read.csv("NAHLN_Data/NAHLN_CattleEcoli_PPY1-PPY5_10-12-2023.csv")
# incidence_matrix <- NAHLN_process_genotype_data(gene_df = gene_df, Year_lookup_df = Year_lookup_df, ID_lookup_df = ID_lookup_df, data_source = "NAHLN", output_level = "gene")
# print(incidence_matrix)