library(stringr)
library(rlist)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(readxl)
library(cowplot)
library(arules)
library(reshape2)
library(gridExtra)
library(grid)
library(lubridate)
library(tidyr)
library(purrr)
library(dplyr)
# library(multcomp)

#function to get a data frame that has genes in one column and corresponding gene classes in another column
get_gene_classes <- function() {

cecal_df <- as.data.frame(read.csv("cecal/cecal_corrected_genotype.csv"))

retail_meats_df <- as.data.frame(read.csv("Retail_Meats/Retail_Meats_corrected_genotype.csv"))

df_list <- list(cecal_df, retail_meats_df)

#make a list to hold each class df
class_df_list <- list()

for (df in df_list) {
df <- df[, which(colnames(df) == "Aminoglycoside_Resist_Genes"):which(colnames(df) == "Other_Classes_Resist_Genes")]

class_list <- list()

for (col in colnames(df)) {
col_string <- c(unlist(strsplit(ifelse(is.na(df[, col]), "", as.character(df[, col])), " ")))
col_string <- unique(col_string[! col_string%in% c(NA, "")]) #get unique genes
col_df <- data.frame(Class = c(replicate(length(col_string), col)), 
            Genes = col_string) #create df



class_list <- list.append(class_list, col_df)
}
class_df <- bind_rows(class_list)
class_df_list <- list.append(class_df_list, class_df)
}
#combine all dfs in list
combined_class_df <- bind_rows(class_df_list)
#get rid of duplicate rows
combined_class_df <- combined_class_df[!duplicated(combined_class_df), ]


#check Other_Classes_Resist_Genes rows, and only keep Genes that do not occur elsewhere in the df
#combined_class_df <- combined_class_df[!(combined_class_df$Class == "Other_Classes_Resist_Genes" & duplicated(combined_class_df$Genes)), ]




write.csv(combined_class_df, file = "gene_class_mappings.csv", row.names = FALSE)
return(combined_class_df)
}











#function to get data into a format that can be used by the apriori algorithm
get_transactions <- function(df, resistance_indicator, year, agg = FALSE) {


    if (resistance_indicator == "phenotype") {
            df <- as.data.frame(df)

            df <- df[df$Year == year, ] #select which year to use



            df <- subset(x = df, select = -c(ID, Year))
            #df <- df[,colSums(is.na(df))<nrow(df)]
            #df <- df[ , ] == 1
            df <- mutate_if(df, is.numeric, as.logical)

            #write.csv(x = df, file = str_glue("incidence_matrices/phenotype/{year}_transactions.csv"), row.names = FALSE)

            trans <- transactions(df) #make a transaction matrix from wide format
            if (agg == TRUE) {
            #add classes to the itemInfo so they can be used for aggregation
            class_mappings <- read_excel("EcoliBreakPoints.xlsx")
            class_list <- list()
            for (label in itemLabels(trans)) {
                for (row in 1:nrow(class_mappings)){

                    if (class_mappings[row, "abbreviation"] == label) {
                        class_list <- list.append(class_list, class_mappings[row, "class"])
                    }
                }
            }
            itemInfo(trans) <- data.frame(labels = itemLabels(trans), class = unlist(class_list))

            
                trans <- aggregate(trans, by = "class")
            }
            
            return(trans)
    }
    if (resistance_indicator == "genotype") {
        
        #df <- dup_genecols_for_cross_res(df)

        df <- df[df$Year == year, ] #select which year to use

        df <- subset(x = df, select = -c(ID, Year))
        df <- df[,colSums(is.na(df))<nrow(df)]
        df <- df[,] == 1

        trans <- transactions(df) #make a transaction matrix from wide format

        if (agg == TRUE) {
        #add classes to the itemInfo so they can be used for aggregation
        class_mappings <- read.csv("gene_class_mappings.csv")

        class_mappings$Gene_family <- gsub(pattern = "[^A-Za-z0-9]", replacement = ".", class_mappings$Genes)
        class_mappings <- class_mappings[!duplicated(class_mappings$Gene_family), ]

        class_list <- list()
        for (label in itemLabels(trans)) {
            for (row in 1:nrow(class_mappings)){

                if (class_mappings[row, "Gene_family"] == label) {
                class_list <- list.append(class_list, class_mappings[row, "corrected_class"])

                    }
                }
            }


            itemInfo(trans) <- data.frame(labels = itemLabels(trans), class = unlist(class_list))

          
                trans <- aggregate(trans, by = "class")
            }



        return(trans)
    }
    else {stop("ERROR: resistance_indicator must be either phenotype or genotype")}

}









#function that gets rules and fequent itemsets by running the apriori algorithm
get_rules_or_itemsets <- function(df, resistance_indicator, year, target, agg = FALSE) {


    if (resistance_indicator == "phenotype") {

        trans <- get_transactions(df = df, resistance_indicator = resistance_indicator,
                                year = year, agg = agg)

        year_df <- df[df$Year == year, ]
        minsupp <- (1 / length(c(unique(year_df$ID)))) #get minimum support


        #get rules/frequent itemsets
        rules_or_itemsets <- apriori(trans,
                    parameter = list(support = minsupp, confidence = 0.0, target = target),
                    control = list(verbose = FALSE))

        return(rules_or_itemsets) #return the rules/frequent itemsets

    }
    if (resistance_indicator == "genotype") {

        trans <- get_transactions(df = df, resistance_indicator = resistance_indicator, year = year, agg = agg)
        year_df <- df[df$Year == year, ]
        minsupp <- (1 / length(c(unique(year_df$ID)))) #get minimum support
        
       
        #get rules
        rules_or_itemsets <- apriori(trans,
                    parameter = list(support = minsupp, confidence = 0.0, target = target),
                    control = list(verbose = FALSE))

        return(rules_or_itemsets) #return the rules or frequent itemsets

    }
    else {stop("ERROR: resistance_indicator must be either phenotype or genotype")}



}




###################################################################################################################
###################################################################################################################
#################################  FUNCTIONS FOR PREPROCESSING  ###################################################
###################################################################################################################
###################################################################################################################

#function to add breakpoints if needed
add_break_points <- function(df, data_source) {
print("STARTING BP ADDITION")
break_points <- as.data.frame(read_excel("EcoliBreakPoints.xlsx"))

# Ensure df is a data frame
df <- as.data.frame(df)

df <- df %>% mutate(across(where(is.character), str_trim)) # get rid of extra white space
df$Drugcode <- toupper(df$Drugcode) # make everything uppercase
df$susceptible <- NA


for (row in 1:nrow(df)) {

    for (bp in 1:nrow(break_points)) {
    if (df[row, "Drugcode"] == break_points[bp, "abbreviation"]) {
        df[row, "resistant"] <- break_points[bp, "resistant_value"]

    }
    }
}

print("ALL DONE ADDING BP")
return(df)

}



#function to add resistance status based on susceptibility breakpoints
add_resistance_status <- function(df, data_source) {
    
print("STARTING RES STATUS ADDITION")

    #fix the min and max where NULL
    #df$Maximum <- df$Maximum %>% replace(.=="NULL", 99999)
    #df$Maximum <- as.numeric(df$Maximum)
    #df$Minimum <- df$Minimum %>% replace(.=="NULL", 0)
    #df$Minimum <- as.numeric(df$Minimum)

    #df <- df %>% mutate(across(where(is.character), str_trim)) # get rid of extra white space
    #df$Drugcode <- toupper(df$Drugcode) # make everything uppercase


    #df <- df[is.finite(df$susceptible), ] #drop NA row-wise

    #loop through each row and determine resistance status
    for (row in 1:nrow(df)) {
    if (df[row, "Minimum"] >= df[row, "resistant"]) {
        df[row, "isResistant"] <- 1
        df[row, "Drugcode"] <- df[row, "Drugcode"]
        #df[row, "Class"] <- df[row, "Class"]



    }
    else {
        df[row, "isResistant"] <- 0
    }

    if(df[row, "Drugcode"] == "FIS"){
        
        if (df[row, "Minimum"] > df[row, "resistant"]) {
        df[row, "isResistant"] <- 1
        df[row, "Drugcode"] <- df[row, "Drugcode"]
        #df[row, "Class"] <- df[row, "Class"]



    }
    else {
        df[row, "isResistant"] <- 0
    }

    }
}







df <- pivot_wider(data = df, id_cols = c(ID, Year), names_from = Drugcode, values_from = isResistant)



write.csv(df, file = str_glue("dataset_specific_outputs/{data_source}/NEW_{data_source}_wide_resStatus_phenotype.csv"), row.names = FALSE)
print("DONE ADDING RES STATUS")
return(df) #return the df with a  isResistant column filled, in wide format

}






#function to check and see if there are any uninterpretable MIC values in the dataset
get_uninterpretable <- function(path) {

    df <- add_resistance_status(path) #read in the data
    uninterpretable_df <- df[(df$Minimum < df$susceptible & df$Maximum > df$SusVal),] #get df with only uninterpretable values

    print(nrow(uninterpretable_df))#print number of uninterpretable values
    write.csv(uninterpretable_df, "uninterpretable_items.csv", row.names = FALSE)#save for further examination

}








prevalence_descriptives <- function(df, data_source, resistance_indicator, class_level = FALSE) {
    # If this is genotype data, create a mapping of original to sanitized names
    original_names_lookup <- read.csv("gene_class_mappings.csv")
    original_names_lookup <- original_names_lookup$Genes

   original_names <- NULL
    if (resistance_indicator == "genotype") {
        # Create sanitized versions of the lookup names to match with df column names
        sanitized_lookup <- make.names(original_names_lookup)
        # Create mapping from sanitized to original names
        name_mapping <- setNames(original_names_lookup, sanitized_lookup)
    }



    year_list <- list()
    sample_size_list <- list()
    df_list <- list()
    
    # Get all possible years in the range
    all_years <- seq(min(df[, "Year"]), max(df[, "Year"]))
    
    for (year in all_years) {
        year_df <- df[df$Year == year, ]
        
        # If there's no data for this year, add NA values
        if (nrow(year_df) == 0) {
            year_list <- list.append(year_list, year)
            sample_size_list <- list.append(sample_size_list, 0)
            # Create a row of NAs for this year
            means_df <- rep(NA, ncol(df) - 2)  # -2 for ID and Year columns
            names(means_df) <- colnames(df)[!colnames(df) %in% c("ID", "Year")]
        } else {
            year_list <- list.append(year_list, unique(year_df$Year))
            sample_size_list <- list.append(sample_size_list, length(unique(year_df$ID)))
            year_df <- subset(year_df, select = -c(ID, Year))
            means_df <- round(x = colMeans(year_df, na.rm = TRUE), digits = 4)
        }
        
        df_list <- list.append(df_list, means_df)
    }

    prevalence_df <- bind_rows(df_list)
    prevalence_df <- prevalence_df * 100
    
    # Sort the columns alphabetically (using sanitized names for genotype data)
    drug_cols <- sort(colnames(prevalence_df))
    prevalence_df <- prevalence_df[, drug_cols]
    
    # Add back the Year and N columns at the beginning
    prevalence_df <- cbind(
        Year = unlist(year_list),
        N = unlist(sample_size_list),
        prevalence_df
    )

    # Calculate summary statistics
    drug_cols <- colnames(prevalence_df)[!(colnames(prevalence_df) %in% c("Year", "N"))]
    
    # Create two summary rows
    mean_row <- data.frame(
        Year = paste(min(prevalence_df$Year), "-", max(prevalence_df$Year)),
        N = sum(prevalence_df$N),
        matrix(NA, nrow = 1, ncol = length(drug_cols), 
               dimnames = list(NULL, drug_cols))
    )
    
    sd_row <- data.frame(
        Year = "SD",
        N = NA,
        matrix(NA, nrow = 1, ncol = length(drug_cols), 
               dimnames = list(NULL, drug_cols))
    )
    
    # Calculate mean and SD for each drug column
    for (col in drug_cols) {
        mean_val <- round(mean(prevalence_df[[col]], na.rm = TRUE), 2)
        sd_val <- round(sd(prevalence_df[[col]], na.rm = TRUE), 2)
        mean_row[[col]] <- sprintf("%.2f", mean_val)
        sd_row[[col]] <- sprintf("%.2f", sd_val)
    }
    
    # Add both summary rows to prevalence_df
    prevalence_df <- rbind(prevalence_df, mean_row, sd_row)

    # If this is genotype data, restore the original column names
    if (resistance_indicator == "genotype") {
        # Get the current column names (excluding Year and N)
        current_cols <- colnames(prevalence_df)[!(colnames(prevalence_df) %in% c("Year", "N"))]
        
        # Map sanitized names back to original names
        new_cols <- c("Year", "N")
        for (col in current_cols) {
            if (col %in% names(name_mapping)) {
                new_cols <- c(new_cols, name_mapping[col])
            } else {
                new_cols <- c(new_cols, col)  # Keep original if no mapping found
            }
        }
        colnames(prevalence_df) <- new_cols
    }

    if (class_level == FALSE) {
        write.csv(prevalence_df, str_glue("dataset_specific_outputs/{data_source}/prevalence_descriptives/{data_source}_{resistance_indicator}_prevalence_descriptives.csv"), row.names = FALSE)
    } else if (class_level == TRUE) {
        write.csv(prevalence_df, str_glue("dataset_specific_outputs/{data_source}/prevalence_descriptives/{data_source}_{resistance_indicator}_class_level_prevalence_descriptives.csv"), row.names = FALSE)
    }
    
    return(prevalence_df)
}















class_level_prevalence_descriptives <- function(df, data_source, resistance_indicator) {

    # CONVERT TO CLASS LEVEL
    if (resistance_indicator == "phenotype") {

        class_mappings <- read_excel("EcoliBreakPoints.xlsx")
    
        df <- pivot_longer(df, cols = -c(ID, Year), names_to = "Drugcode",  values_to = c("isResistant"))

        code_vec <- c(class_mappings$abbreviation)
        class_vec <- c(class_mappings$class)

        for (row in 1:nrow(df)) {
            for (code in 1:length(code_vec)) {
                if (code_vec[code] == df[row, "Drugcode"]) {
                    df[row, "Drugcode"] <-  class_vec[code]
                }
            }
        }

        # Take care of the duplicates problem (results from aggregating multiple drugs into the same class)
        df <- arrange(df, Year, ID, Drugcode, -isResistant)
        df <- df[!base::duplicated(df[, colnames(df)[colnames(df) != 'isResistant']]),]

        df <- pivot_wider(data = df, id_cols = c(ID, Year), names_from = Drugcode, values_from = isResistant)
    }

    if (resistance_indicator == "genotype") {
        print("original colnames:")
        print(colnames(df))
        genotype_class_mappings <- read.csv("gene_class_mappings.csv")
        genotype_class_mappings$Genes <- gsub(pattern = "[^A-Za-z0-9]", replacement = ".", genotype_class_mappings$Genes)
        
        # Ensure corrected_class is a character vector
        genotype_class_mappings$corrected_class <- as.character(genotype_class_mappings$corrected_class)

        genotype_class_names <- setNames(make.unique(genotype_class_mappings$corrected_class), genotype_class_mappings$Genes)
        print("Gene class mappings:")
        print(genotype_class_names)
        genotype_duplicated_classes <- names(which(table(genotype_class_mappings$corrected_class) > 1))

        # Function to merge duplicated columns keeping the largest value
        merge_columns <- function(df, cols) {
            col_name <- unique(gsub("\\.\\d+$", "", cols[1]))
            df %>% 
                mutate(!!col_name := pmax(!!!rlang::syms(cols), na.rm = TRUE)) %>%
                select(-all_of(cols), !!col_name)
        }

        df <- df %>%
            # First, rename the columns
            rename_with(~ genotype_class_names[.x], .cols = intersect(names(.), names(genotype_class_names))) %>%
            {   
                df <- .
                # Then, merge the duplicated columns
                for (class in genotype_duplicated_classes) {
                    cols_to_merge <- names(df)[grepl(paste0("^", class, "(\\.|$)"), names(df))]
                    print(paste("For class:", class))
                    print(paste("Columns to merge:", paste(cols_to_merge, collapse=", ")))
                    if (length(cols_to_merge) > 1) {
                        df <- merge_columns(df, cols_to_merge)
                    }
                }
                df
            }

        # Remove trailing dot and number from column names in genotype_df
        colnames(df) <- gsub("\\.\\d+$", "", colnames(df))

        # Check for unmapped genes
        unmapped_genes <- setdiff(names(df), c("ID", "Year", unique(genotype_class_mappings$corrected_class)))
        if (length(unmapped_genes) > 0) {
            warning("The following genes were not mapped to any class: ", paste(unmapped_genes, collapse = ", "))
        }

        # Debugging: Print the final column names
        print("Final column names:")
        print(colnames(df))
    }

    # GET THE PREVALENCE DESCRIPTIVES
    year_list <- list()
    sample_size_list <- list()
    df_list <- list()
    for (year in c(min(df[, "Year"]):(max(df[, "Year"])))){

        year_df <- df[df$Year == year,]
        year_list <- list.append(year_list, unique(year_df$Year))
        sample_size_list <- list.append(sample_size_list, length(unique(year_df$ID)))
        year_df <- subset(year_df, select = -c(ID, Year))
        
        means_df <- round(x = colMeans(year_df, na.rm = TRUE), digits = 4)

        df_list <- list.append(df_list, means_df)
    }

    prevalence_df <- bind_rows(df_list)
    prevalence_df <- prevalence_df * 100
    prevalence_df <- cbind(N = unlist(sample_size_list), prevalence_df)
    prevalence_df <- cbind(Year = unlist(year_list), prevalence_df)

    write.csv(prevalence_df, str_glue("dataset_specific_outputs/{data_source}/prevalence_descriptives/{data_source}_{resistance_indicator}_class_level_prevalence_descriptives.csv"), row.names = FALSE)
}









#function to get the number of isolates that have a gene that is not represented in the phenotype data at the class level
isolate_level_gene_phenotype_mismatch <- function(genotype_df, phenotype_df, data_source) {




    phenotype_class_mappings <- read_excel("EcoliBreakPoints.xlsx")
    genotype_class_mappings <- read.csv("gene_class_mappings.csv")
    genotype_class_mappings$Genes <- gsub(pattern = "[^A-Za-z0-9]", replacement = ".", genotype_class_mappings$Genes)


    # Create a named vector for renaming, ensuring uniqueness
    phenotype_class_names <- setNames(make.unique(phenotype_class_mappings$class), phenotype_class_mappings$abbreviation)
    genotype_class_names <- setNames(make.unique(genotype_class_mappings$corrected_class), genotype_class_mappings$Genes)

    # Identify duplicated class names (based on original names, before make.unique)
    phenotype_duplicated_classes <- names(which(table(phenotype_class_mappings$class) > 1))
    genotype_duplicated_classes <- names(which(table(genotype_class_mappings$corrected_class) > 1))

    # Function to merge duplicated columns keeping the largest value
    merge_columns <- function(df, cols) {
        col_name <- unique(gsub("\\.\\d+$", "", cols[1]))
        df %>% 
            mutate(!!sym(col_name) := pmax(!!!syms(cols), na.rm = TRUE)) %>%
            select(-all_of(cols)) %>%
            select(everything(), !!sym(col_name))  # Move the new column to the end
    }

    # Rename and merge columns in phenotype_df
    phenotype_df <- phenotype_df %>%
        # First, rename the columns
        rename_with(~ phenotype_class_names[.x], .cols = intersect(names(.), names(phenotype_class_names))) %>%
        {
            df <- .
            # Then, merge the duplicated columns
            for (class in phenotype_duplicated_classes) {
                cols_to_merge <- names(df)[grepl(paste0("^", class, "(\\.|$)"), names(df))]
                if (length(cols_to_merge) > 1) {
                    df <- merge_columns(df, cols_to_merge)
                }
            }
            df
        }
    # Rename and merge columns in genotype_df
    genotype_df <- genotype_df %>%
        # First, rename the columns
        rename_with(~ genotype_class_names[.x], .cols = intersect(names(.), names(genotype_class_names))) %>%
        {   
            df <- .
            # Then, merge the duplicated columns
            for (class in genotype_duplicated_classes) {
                cols_to_merge <- names(df)[grepl(paste0("^", class, "(\\.|$)"), names(df))]
                
                if (length(cols_to_merge) > 1) {
                    df <- merge_columns(df, cols_to_merge)
                }
            }
            df
        }
   # Remove trailing dot and number from column names in genotype_df
   colnames(genotype_df) <- gsub("\\.\\d+$", "", colnames(genotype_df))

    #remove "-ECEC" from the end of the ID in genotype_df
    genotype_df$ID <- gsub("-ECEC", "", genotype_df$ID)

    #also remove "EC" from the end of the ID in genotype_df
    genotype_df$ID <- gsub("EC", "", genotype_df$ID)


    #Include only rows where the ID matches in both dataframes
    phenotype_df <- phenotype_df[phenotype_df$ID %in% genotype_df$ID, ]
    genotype_df <- genotype_df[genotype_df$ID %in% phenotype_df$ID, ]
   
    

    #drop ID (i.e., keep only the Year andclass level antimicrobial resistance status columns)
    #phenotype_df <- subset(phenotype_df, select = -c(ID))
    #genotype_df <- subset(genotype_df, select = -c(ID))

  
    
    #make Year a factor
    phenotype_df$Year <- as.factor(phenotype_df$Year)
    genotype_df$Year <- as.factor(genotype_df$Year)

    #makes sure the columns are in the same order for both dataframes
    #phenotype_df <- phenotype_df[, colnames(genotype_df)]
    phenotype_df <- phenotype_df[, sort(names(phenotype_df))]
    genotype_df <- genotype_df[, sort(names(genotype_df))]





    # Create a list to store mismatches for each year
    mismatch_list <- list()

    for (year in unique(phenotype_df$Year)) {
        year_phenotype <- phenotype_df[phenotype_df$Year == year, ]
        year_genotype <- genotype_df[genotype_df$Year == year, ]
        
        # Count mismatches, including NAs
        mismatches <- colSums(
            (is.na(year_phenotype) & !is.na(year_genotype)) |
            (!is.na(year_phenotype) & is.na(year_genotype)) |
            (year_phenotype != year_genotype & !is.na(year_phenotype) & !is.na(year_genotype))
        )
        mismatch_list[[as.character(year)]] <- mismatches
    }

    # Convert the list to a data frame
    mismatch_df <- do.call(rbind, mismatch_list)

    # Add row names (years) and column names
    rownames(mismatch_df) <- unique(phenotype_df$Year)
    colnames(mismatch_df) <- colnames(phenotype_df)

    # Remove the 'Year' and 'ID' columns if they exist
    mismatch_df <- mismatch_df[, !colnames(mismatch_df) %in% c("Year", "ID")]



    # Convert mismatch_df to a data frame with row names as a column
    mismatch_table <- data.frame(Year = rownames(mismatch_df), mismatch_df, check.names = FALSE)

# Function to create a themed table grob
create_table_grob <- function(data, title) {
  table_theme <- ttheme_minimal(
    core = list(fg_params = list(fontface = c(rep("plain", ncol(data) - 1), "bold")),
                bg_params = list(fill = c(rep(c("white", "grey95"), length.out = nrow(data))))),
    colhead = list(fg_params = list(fontface = "bold"),
                   bg_params = list(fill = "grey80")),
    rowhead = list(fg_params = list(fontface = "bold"))
  )
  
  grob <- tableGrob(data, rows = NULL, theme = table_theme)
  
  # Add title
  title_grob <- textGrob(title, gp = gpar(fontface = "bold", fontsize = 20))
  padding <- unit(1.5, "line")
  
  grob <- gtable::gtable_add_rows(grob, heights = grobHeight(title_grob) + padding, pos = 0)
  grob <- gtable::gtable_add_grob(grob, title_grob, t = 1, l = 1, r = ncol(grob))
  
  return(grob)
}

# Create the table grob
table_grob <- create_table_grob(mismatch_table, "Frequency of Isolate Level Mismatches (Genotype vs Phenotype)")

# Save the table as a PNG file
png("mismatch_table.png", width = 10, height = 3, units = "in", res = 300)
grid.draw(table_grob)
dev.off()

    # Return the mismatch dataframe
    return(mismatch_df)
}













compare_dataframes <- function(df1, df2) {
  # Check if dimensions are the same
  if (!all(dim(df1) == dim(df2))) {
    print("FALSE")
  }
  
  # Sort column names
  df1 <- df1[, sort(names(df1))]
  df2 <- df2[, sort(names(df2))]
  print("df1 col names")
  print(colnames(df1))
  print("df2 col names")
  print(colnames(df2))

  # Check if column names are the same
  if (!all(names(df1) == names(df2))) {
    print("FALSE")
  }
  
  # Compare values
are_equivalent <- all(df1 == df2, na.rm = TRUE)
print(paste("Dataframes are equivalent:", are_equivalent))

# If not equivalent, find differences
if (!are_equivalent) {
  differences <- which(df1 != df2, arr.ind = TRUE)
  print("Differences found at:")
  print(differences)
  
  # Print a few example differences
  for (i in 1:min(20, nrow(differences))) {
    row <- differences[i, 1]
    col <- differences[i, 2]
    cat(sprintf("Row %d, Column '%s': df1 = %s, df2 = %s\n", 
                row, names(df1)[col], df1[row, col], df2[row, col]))
  }
}

}












#function to match the gene family from refgenes with genes from NARMS data
match_to_refgenes <- function(AMRFinder_df, refgenes_df) {
 # Initialize data frame with the correct number of rows
 narms_df <- data.frame(
   Genes = unique(AMRFinder_df$Element.symbol),
   stringsAsFactors = FALSE
 )
  # Get rid of rows where the Gene.family column is duplicated
  refgenes_df <- refgenes_df[!duplicated(refgenes_df$Gene.family), ]
  
 find_closest_match <- function(gene, ref_genes) {
  # replace ` with ' in gene
  clean_gene <- gsub("`", "'", gene)
  ref_gene <- gsub("`", "'", ref_genes)

  # First, try to find an exact match (case-insensitive)
  exact_match <- which(tolower(ref_genes) == tolower(clean_gene))
  if (length(exact_match) > 0) {
    return(list(match = ref_genes[exact_match[1]], distance = 0))
  }
  
  # Next, try to find a match where one string starts with the other
  starts_with_match <- which(startsWith(tolower(ref_genes), tolower(clean_gene)) | 
                             startsWith(tolower(clean_gene), tolower(ref_genes)))
  if (length(starts_with_match) > 0) {
    return(list(match = ref_genes[starts_with_match[1]], distance = 0.1))  # Use 0.1 to indicate a very close match
  }
  
  # If no exact or starts-with match is found, use Levenshtein distance
  distances <- stringdist::stringdist(clean_gene, ref_genes, method = "lv")
  closest <- which.min(distances)
  
  # Only return a match if the Levenshtein distance is less than half the length of the shorter string
  min_length <- min(nchar(clean_gene), nchar(ref_genes[closest]))
  if (distances[closest] <= min_length / 2) {
    return(list(match = ref_genes[closest], distance = distances[closest]))
  } else {
    return(list(match = NA, distance = Inf))  # No good match found
  }
}

  
  # Create a new list to store renamed columns
  new_gene_names <- list()
  
  # Create a list to store non-exact matches
  non_exact_matches <- list()
  
  # Iterate through Genes of narms_df
  #check against refgenes_df$#Allele first and then check against refgenes_df$Gene.family
# Iterate through Genes of narms_df
# Check against refgenes_df$Allele first and then check against refgenes_df$Gene.family
for (row in 1:nrow(narms_df)) {
    match_result <- find_closest_match(narms_df[row, "Genes"], refgenes_df$Gene.family)
    closest_match <- match_result$match
    narms_df[row, "Gene_family"] <- closest_match
    narms_df[row, "match_distance"] <- match_result$distance
    
    # Find the matching row in refgenes_df
    matching_row <- which(refgenes_df$Gene.family == closest_match)
    
    # Only assign values if a matching row is found
    if (length(matching_row) > 0) {
        narms_df[row, "corrected_class"] <- refgenes_df$Class[matching_row[1]]
        narms_df[row, "sub_class"] <- refgenes_df$Subclass[matching_row[1]]
    } else {
        narms_df[row, "corrected_class"] <- NA
        narms_df[row, "sub_class"] <- NA
    }
    
    if (match_result$distance > 0) {
        non_exact_matches[[narms_df[row, "Genes"]]] <- closest_match
    }
}

  print("MADE IT HERE")


  cat("Exact matches:\n")
  for (original in names(new_gene_names)) {
    cat(sprintf("%s -> %s\n", original, new_gene_names[[original]]))
  }
  
  # Print non-exact matches
  cat("Non-exact matches:\n")
  for (original in names(non_exact_matches)) {
    cat(sprintf("%s -> %s\n", original, non_exact_matches[[original]]))
  }
  write.csv(narms_df, file = "gene_class_mappings.csv", row.names = FALSE)
  #return(result_df)
}
# AMRFinder_df <- read.csv("Gene_Data/amrfinder_results-AMR-CORE.point-muts-included.csv")
# refgenes_df <- read.csv("refgenes.csv")

# match_to_refgenes(AMRFinder_df, refgenes_df)















convert_gene_to_gene_family <- function(df, genotype_class_mappings, data_source) {
  # Create helper function to standardize names for matching
  standardize_name <- function(name) {
    # Replace any character that isn't a letter or number with a period
    # Also converts to lowercase for case-insensitive matching
    gsub("[^A-Za-z0-9]", ".", tolower(name))
  }
  
  # Standardize names in both dataframes for matching
  std_df_cols <- standardize_name(names(df))                # Standardize column names from input df
  std_genes <- standardize_name(genotype_class_mappings$Genes)  # Standardize genes from mapping file
  
  # Create mapping dictionary from standardized gene names to gene families
  # setNames creates a named vector where:
  # - The values are the gene families
  # - The names are the standardized gene names
  gene_to_family <- setNames(
    genotype_class_mappings$Gene_family,
    std_genes
  )
  
  # Initialize result dataframe with just ID and Year columns
  result <- df[, c("ID", "Year")]
  
  # Process each gene family one at a time
  for (family in unique(genotype_class_mappings$Gene_family)) {
    # Find all standardized gene names that map to this family
    family_genes <- names(gene_to_family)[gene_to_family == family]
    
    # Find which columns in the input df correspond to these genes
    matching_cols <- which(std_df_cols %in% family_genes)
    
    if (length(matching_cols) > 0) {
      # For this gene family:
      # 1. Select all columns that correspond to genes in this family
      # 2. Calculate rowSums to see if ANY gene is present (sum > 0)
      # 3. Convert to numeric to get 1/0 instead of TRUE/FALSE
      result[[family]] <- as.numeric(rowSums(df[, matching_cols, drop = FALSE]) > 0)
    }
  }
  
  # Reorder columns to ensure ID and Year are first, followed by the rest in alphabetical order
  result <- result[, c("ID", "Year", sort(names(result)[-(1:2)]))]

  # Save the result to a CSV file
  write.csv(result, file = str_glue("dataset_specific_outputs/{data_source}/{data_source}_wide_family_level_genotype.csv"), row.names = FALSE)
  
  return(result)
}
















#function to convert the phenotype data to the class level
convert_phenotype_to_class <- function(df, phenotype_class_mappings, data_source) {
  
  # Function to merge case-insensitive duplicate columns
  merge_case_insensitive_columns <- function(df) {
  # Get all column names
  col_names <- names(df)
  
  # Create a list to store columns to be merged
  to_merge <- list()
  
  # Identify columns to be merged
  for (name in col_names) {
    lower_name <- tolower(name)
    if (!(lower_name %in% names(to_merge))) {
      to_merge[[lower_name]] <- col_names[tolower(col_names) == lower_name]
    }
  }
  
  # Merge identified columns
  for (lower_name in names(to_merge)) {
    cols <- to_merge[[lower_name]]
    if (length(cols) > 1) {
      # Keep the capitalization that appears first
      kept_name <- cols[1]
      df[[kept_name]] <- apply(df[, cols, drop = FALSE], 1, max, na.rm = TRUE)
      df <- df[, !(names(df) %in% cols[-1])]
    }
  }
  
  return(df)
}
  
  # Merge case-insensitive duplicate columns
  df <- merge_case_insensitive_columns(df)
print("original column names:")
print(colnames(df))

 #phenotype_class_mappings <- read.csv("gene_class_mappings.csv")

        phenotype_class_names <- setNames(make.unique(phenotype_class_mappings$class), phenotype_class_mappings$abbreviation)
        print("Phenotype class mappings:")
        print(phenotype_class_names)
        phenotype_duplicated_classes <- names(which(table(phenotype_class_mappings$class) > 1))

        # Function to merge duplicated columns keeping the largest value
        merge_columns <- function(df, cols) {
            col_name <- unique(gsub("\\.\\d+$", "", cols[1]))
            df %>% 
                mutate(!!col_name := pmax(!!!rlang::syms(cols), na.rm = TRUE)) %>%
                select(-all_of(cols), !!col_name)
        }

      #merge any columns in df that are currently duplicated
      for (col in colnames(df)) {
                    cols_to_merge <- names(df)[grepl(paste0("^", col, "(\\.|$)"), names(df))]
                    
                    if (length(cols_to_merge) > 1) {
                        df <- merge_columns(df, cols_to_merge)
                    }
                }


        df <- df %>%
            # First, rename the columns
            rename_with(~ phenotype_class_names[.x], .cols = intersect(names(.), names(phenotype_class_names))) %>%
            {   
                df <- .
                # Then, merge the duplicated columns
                for (class in phenotype_duplicated_classes) {
                    cols_to_merge <- names(df)[grepl(paste0("^", class, "(\\.|$)"), names(df))]
                    
                    if (length(cols_to_merge) > 1) {
                        df <- merge_columns(df, cols_to_merge)
                    }
                }
                df
            }

        # Remove trailing dot and number from column names in genotype_df
        colnames(df) <- gsub("\\.\\d+$", "", colnames(df))

        # Check for unmapped genes
        unmapped_drugs <- setdiff(names(df), c("ID", "Year", unique(phenotype_class_mappings$class)))
        if (length(unmapped_drugs) > 0) {
            warning("The following drugs were not mapped to any class: ", paste(unmapped_drugs, collapse = ", "))
        }

        # Function to split columns with "/"
        split_slash_columns <- function(df) {
            slash_cols <- grep("/", names(df), value = TRUE)
            for (col in slash_cols) {
                parts <- strsplit(col, "/")[[1]]
                for (part in parts) {
                    if (!part %in% names(df)) {
                        df[[part]] <- df[[col]]
                    } else {
                        df[[part]] <- pmax(df[[part]], df[[col]], na.rm = TRUE)
                    }
                }
                df[[col]] <- NULL
            }
            return(df)
        }

        # Apply the split_slash_columns function
        df <- split_slash_columns(df)

        # Debugging: Print the final column names
        print("Final column names:")
        print(colnames(df))
        write.csv(df, file = str_glue("dataset_specific_outputs/{data_source}/NEW_{data_source}_wide_class_level_phenotype.csv"), row.names = FALSE)
        
        return(df)
    }





