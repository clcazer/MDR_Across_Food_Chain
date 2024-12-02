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

wrangle_NAHLN_phenotype <- function() {
    df <- read_excel("NAHLN_CattleEcoli_PPY1-PPY5_10-12-2023.xlsx")
    break_points_df <- read_excel("EcoliBreakPoints.xlsx")
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




wrangle_NAHLN_genotype <- function() {
    df <- read_excel("Isolate+BioSample.xlsx")
    df_time <- read_excel("NAHLN_CattleEcoli_PPY1-PPY5_10-12-2023.xlsx")
    #replace spaces in colnames with periods
    colnames(df) <- gsub(" ", ".", colnames(df))

    # Convert TestDate to Date type and extract year
    df_time$Year <- lubridate::year(as.Date(df_time$TestDate))
    
    # Print sample of data for debugging
    cat("Sample of df$Isolate.identifiers:\n")
    print(head(df$Isolate.identifiers))
    cat("\nSample of df_time$UniqueID:\n")
    print(head(df_time$UniqueID))
    
    # Create a mapping of IDs to years and UniqueIDs
    id_to_year <- numeric(length = nrow(df))
    id_to_uniqueid <- character(length = nrow(df))
    
    for (i in 1:nrow(df)) {
        # Find which UniqueID from df_time appears in this Isolate.identifiers
        matching_row <- which(sapply(df_time$UniqueID, function(id) {
            # Ensure both strings are character type and use fixed=TRUE for exact matching
            grepl(as.character(id), as.character(df$Isolate.identifiers[i]), fixed=TRUE)
        }))
        
        if (length(matching_row) > 0) {
            id_to_year[i] <- df_time$Year[matching_row[1]]
            id_to_uniqueid[i] <- df_time$UniqueID[matching_row[1]]
            cat(sprintf("Matched: %s with %s\n", df$Isolate.identifiers[i], df_time$UniqueID[matching_row[1]]))
        } else {
            warning(paste("No matching UniqueID found for Isolate.identifiers:", df$Isolate.identifiers[i]))
        }
    }
    
    # Add Year and UniqueID to df
    df$Year <- id_to_year
    df$UniqueID <- id_to_uniqueid
    
    # Remove rows where Year is 0 instead of reassigning them
    df <- df[df$Year != 0, ]

    # Print summary of year assignments
    cat("\nYear distribution:\n")
    print(table(df$Year))
    
    #replace spaces in colnames with periods
    colnames(df) <- gsub(" ", ".", colnames(df))
    
    df <- df[, c("Isolate.identifiers", "UniqueID", "Year", "AMR.genotypes")]
    
    # Create new column with clean gene names, only keeping COMPLETE genes
    df$clean_genotypes <- sapply(df$AMR.genotypes, function(x) {
        # Split by comma
        genes <- unlist(strsplit(x, ","))
        # Keep only genes with COMPLETE status
        genes <- genes[grepl("=COMPLETE$", genes)]
        # Clean remaining genes and join with spaces
        clean_genes <- gsub("=.*$", "", genes)
        paste(clean_genes, collapse = " ")
    })
    
    # Get all unique genes
    all_genes <- unique(unlist(strsplit(df$clean_genotypes, " ")))
    
    # Create incidence matrix
    incidence_matrix <- matrix(0, 
                             nrow = nrow(df), 
                             ncol = length(all_genes),
                             dimnames = list(df$Isolate.identifiers, all_genes))
    
    # Fill matrix with 1s where genes are present
    for(i in 1:nrow(df)) {
        sample_genes <- unlist(strsplit(df$clean_genotypes[i], " "))
        incidence_matrix[i, sample_genes] <- 1
    }
    
    # Convert to dataframe and add year
    incidence_df <- as.data.frame(incidence_matrix)
    incidence_df$Year <- df$Year
    incidence_df$ID <- df$Isolate.identifiers
    incidence_df$UniqueID <- df$UniqueID
    
    #Reorder columns to put ID first and Year second
    incidence_df <- incidence_df[, c("UniqueID", "ID", "Year", setdiff(names(incidence_df), c("UniqueID","ID", "Year")))]
    
    # Write to CSV
    write.csv(incidence_df, "NAHLN/NAHLN_genotype_wide.csv", row.names = FALSE)
    
    return(incidence_df)
}




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


# df <- read.csv("NAHLN/NAHLN_wide_family_level_genotype.csv")

# top_genes_df <- read.csv("NAHLN_topGenes.csv")


