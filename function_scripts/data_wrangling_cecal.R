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













cecal_MIC_to_interval <- function(df, data_source, limit_date_range = FALSE) {

    #only get the data we want to analyze
    df <- df[df$GROWTH == "YES", ]
    df <- df[df$GENUS == "Escherichia", ]
    df <- df[df$SPECIES == "coli", ]
    df <- df[df$HOST_SPECIES == "Cattle", ]
    df <- df[df$SOURCE == "Cecal (Beef)", ]

   
    if (limit_date_range == TRUE) {
    #only get 2017 onward
    df <- df[df$Year >= 2017, ]
    }
    print(colnames(df))


    #INCLUDE STR if TRUE
    if (limit_date_range == TRUE) {
    #split the df into a sign df and an MIC df (will make pivoting easier)
    sign_df <- df[, c("ID", "Year", "SOURCE",
    "AMC.Sign",
    "AMP.Sign",
    "AXO.Sign",
    "AZI.Sign",
    "CHL.Sign",
    "CIP.Sign",
    "COT.Sign",
    "FIS.Sign",
    "FOX.Sign",
    "GEN.Sign",
    "MER.Sign",
    "NAL.Sign",
    "SMX.Sign",
    "TET.Sign",
    "STR.Sign")]




    MIC_df <- df[, c("ID", "Year", "SOURCE",
    "AMC",
    "AMP",
    "AXO",
    "AZI",
    "CHL",
    "CIP",
    "COT",
    "FIS",
    "FOX",
    "GEN",
    "MER",
    "NAL",
    "SMX",
    "TET",
    "STR")]

    } else { #ELSE EXCLUDE STR
    sign_df <- df[, c("ID", "Year", "SOURCE",
    "AMC.Sign",
    "AMP.Sign",
    "AXO.Sign",
    "AZI.Sign",
    "CHL.Sign",
    "CIP.Sign",
    "COT.Sign",
    "FIS.Sign",
    "FOX.Sign",
    "GEN.Sign",
    "MER.Sign",
    "NAL.Sign",
    "SMX.Sign",
    "TET.Sign")]




    MIC_df <- df[, c("ID", "Year", "SOURCE",
    "AMC",
    "AMP",
    "AXO",
    "AZI",
    "CHL",
    "CIP",
    "COT",
    "FIS",
    "FOX",
    "GEN",
    "MER",
    "NAL",
    "SMX",
    "TET")]

    }
    #check to see if MIC_df and sign df are empty
    if (nrow(sign_df) == 0) {
        print("sign_df is empty")
    }
    if (nrow(MIC_df) == 0) {
        print("MIC_df is empty")
    }



    #pivot the sign_df to long format delete the .Sign from the Drugcode column so we are just left with drug code that will match that on the MIC_df
    sign_df <- pivot_longer(sign_df, cols = ends_with("Sign"), names_to = "Drugcode", values_to = c("MIC_SIGN"), values_drop_na = TRUE)# pivot to long format
    sign_df$Drugcode  <- substr(sign_df$Drugcode , 0, 3)


    #pivot the MIC_df to long format
    MIC_df <- pivot_longer(MIC_df, cols = -c(ID, Year, SOURCE), names_to = "Drugcode", values_to = c("MIC"), values_drop_na = TRUE)# pivot to long format



    #now merge both long formatted dataframes
    merged_df <- merge(x = sign_df, y = MIC_df, fill = NA)


    #add Minimum and Maximum columns to the merged_df
    merged_df$Minimum <- NA
    merged_df$Maximum <- NA


    #loop through merged_df and add values to min and max columns based on MIC_SIGN and MIC column values
    for (row in 1:nrow(merged_df)) {

        if (merged_df[row, "MIC_SIGN"] == "="){

            merged_df[row, "Maximum"] <- merged_df[row, "MIC"]
            merged_df[row, "Minimum"] <- merged_df[row, "MIC"]
        }
        else if (merged_df[row, "MIC_SIGN"] == "<") {

            merged_df[row, "Maximum"] <- as.numeric(merged_df[row, "MIC"]) - 0.1
            merged_df[row, "Minimum"] <- 0
        }
        else if (merged_df[row, "MIC_SIGN"] == ">") {

            merged_df[row, "Maximum"] <- 99999
            merged_df[row, "Minimum"] <- as.numeric(merged_df[row, "MIC"]) + 0.1
        }
        else if (merged_df[row, "MIC_SIGN"] == "<=") {

            merged_df[row, "Maximum"] <- merged_df[row, "MIC"]
            merged_df[row, "Minimum"] <- 0
        }
        else if (merged_df[row, "MIC_SIGN"] == ">=") {

            merged_df[row, "Maximum"] <- 99999
            merged_df[row, "Minimum"] <- merged_df[row, "MIC"]
        }

    }


    write.csv(merged_df, file = str_glue("{data_source}/NEW_{data_source}_phenotype_data.csv"), row.names = FALSE)
    print("DONE MAKING THE MIC INTERVALS")
    return(merged_df)

}





cecal_process_genotype_data <- function(gene_df, ID_Year_lookup_df, data_source, output_level = 'gene', limit_date_range = FALSE) {
    gene_df <- gene_df[gene_df$meta == "NARMS_CECAL", ]

    ID_Year_lookup_df <- ID_Year_lookup_df[ID_Year_lookup_df$SOURCE == "Cecal (Beef)", ]

    print(length(unique(gene_df$SAMN)))
    print(length(gene_df$SAMN))
    
    # Create a mapping dataframe with just the columns we need
    id_year_mapping <- ID_Year_lookup_df[, c("BIOSAMPLE_ID", "SRA_ACCESSION_NUMBER", "ID", "Year")]
    
    # Merge the gene_df with the mapping dataframe
    gene_df <- merge(gene_df, 
                    id_year_mapping, 
                    by.x = "SAMN", 
                    by.y = "SRA_ACCESSION_NUMBER", 
                    all.x = TRUE)
    
    # Print the number of rows before and after to verify the merge
    # print(paste("Number of rows with ID/Year:", sum(!is.na(gene_df$ID))))

    
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

    if (limit_date_range == TRUE) {
        gene_df <- gene_df[gene_df$Year >= 2017, ]
    }
    
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


#  gene_df <- read.csv("Gene_Data/amrfinder_results-AMR-CORE.point-muts-included.csv")
#  ID_Year_lookup_df <- read.csv("Slaughterhouse_Data/old/CVM-2020-NARMS-Cecal-Data.csv")
#  incidence_matrix <- cecal_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "cecal", output_level = "gene")
#  print(incidence_matrix)