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

#get_gene_classes()


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
        print(year)
       
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

break_points <- as.data.frame(read_excel("EcoliBreakPoints.xlsx"))


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

#write.csv(df, str_glue("{data_source}/{data_source}_phenotype_withBP_NOTdb.csv"), row.names = FALSE)
print("ALL DONE ADDING BP")
return(df)

}



#function to add resistance status based on susceptibility breakpoints
add_resistance_status <- function(df, data_source) {
    


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



write.csv(df, file = str_glue("{data_source}/{data_source}_wide_resStatus_phenotype.csv"), row.names = FALSE)
return(df) #return the df with a  isResistant column filled, in wide format

}


#function to check and see if there are any uninterpretable MIC values in the dataset
get_uninterpretable <- function(path) {

    df <- add_resistance_status(path) #read in the data
    uninterpretable_df <- df[(df$Minimum < df$susceptible & df$Maximum > df$SusVal),] #get df with only uninterpretable values

    print(nrow(uninterpretable_df))#print number of uninterpretable values
    write.csv(uninterpretable_df, "uninterpretable_items.csv", row.names = FALSE)#save for further examination

}

get_genes_wide <- function(df, data_source) {
    df <- df[, c("ID","Year","SOURCE","CORRECTED_GENOTYPE")] #only want certain columns
    
    print("Before processing:")
    print(str(df$CORRECTED_GENOTYPE))
    print(head(df$CORRECTED_GENOTYPE, 20))
    
    # Split the genes
    len <- sapply(gregexpr("\\S+", df$CORRECTED_GENOTYPE), length)
    df$len <- len
    df_new <- data.frame(cbind(df, str_split_fixed(df$CORRECTED_GENOTYPE, " ", max(len))))
    names(df_new) <- c("ID", "Year", "SOURCE", "CORRECTED_GENOTYPE", "length", paste0("GENE", 1:max(len)))


    df_new <- subset(df_new, select = -c(SOURCE, CORRECTED_GENOTYPE, length)) #delete unneeded 
    df_new <- pivot_longer(df_new, cols = 3:ncol(df_new)) # pivot to long format
    df_new <- df_new[df_new$value != "", ]
    df_new[, "isPresent"] <- 1

    df_new <- pivot_wider(data = df_new, id_cols = c(ID, Year), names_from = value, values_from = isPresent, values_fill = list(isPresent = 0))

    #drop "name" column if it exists
    if("name" %in% colnames(df_new)) {
        df_new <- df_new[, !colnames(df_new) %in% "name"]
    }

    # Remove "NA" column if it exists
    if("NA" %in% colnames(df_new)) {
        df_new <- df_new[, !colnames(df_new) %in% "NA"]
    }

    #drop NULL column if it exists
    if("NULL" %in% colnames(df_new)) {
        df_new <- df_new[, !colnames(df_new) %in% "NULL"]
    }



# Function to clean column names and merge duplicates
clean_and_merge_columns <- function(df) {
  # Clean column names
  colnames(df) <- gsub(" NA$", "", colnames(df))
  
  # Identify duplicate columns
  col_names <- colnames(df)
  duplicates <- col_names[duplicated(col_names) | duplicated(col_names, fromLast = TRUE)]
  
  # Merge duplicate columns
  for (col in unique(duplicates)) {
    cols_to_merge <- which(col_names == col)
    df[[col]] <- do.call(pmax, c(df[cols_to_merge], na.rm = TRUE))
    df <- df[, -cols_to_merge[-1]]  # Remove all but the first occurrence
  }
  
  return(df)
}



df_new <- clean_and_merge_columns(df_new)

#check length of unique IDs for each year
    for (year in unique(df_new$Year)) {
        df <- df_new[df_new$Year == year, ]
        print(length(unique(df$ID)))
    }

    #remove rows with duplicate IDs
    #df_new <- df_new[!duplicated(df_new$ID), ]

    print("Final column names:")
    print(colnames(df_new))
    write.csv(df_new, file = str_glue("{data_source}/{data_source}_wide_corrected_genotype.csv"), row.names = FALSE)
}


#corrected_retail_meats_df <- read.csv("Retail_Meats/Retail_Meats_corrected_genotype.csv")
#get_genes_wide(corrected_retail_meats_df, "Retail_Meats")


#this function is needed because there are sometimes genes in the genotype column that are not in the class labeled columns and vice versa
fix_genotype <- function(df, data_source) {

#replace all spaces, dashes, and periods in column names with underscores
colnames(df) <- gsub(pattern = " ", replacement = "_", x = colnames(df))
colnames(df) <- gsub(pattern = "-", replacement = "_", x = colnames(df))
colnames(df) <- gsub(pattern = "\\.", replacement = "_", x = colnames(df))


# Filter out rows with empty or NA BIOSAMPLE_ID and SRA_ACCESSION_NUMBER only if GENOTYPE is blank
df <- df[!(df$GENOTYPE == "" | is.na(df$GENOTYPE)) | 
          (!is.na(df$BIOSAMPLE_ID) & df$BIOSAMPLE_ID != "") | 
          (!is.na(df$SRA_ACCESSION_NUMBER) & df$SRA_ACCESSION_NUMBER != ""), ]


#only get the data we want to analyze
df <- df[df$GROWTH == "YES", ]
df <- df[df$GENUS == "EC", ]
df <- df[df$SPECIES == "coli", ]
df <- df[df$HOST_SPECIES == "Cattle", ]
df <- df[df$Year != 2020, ]
df <- df[df$Year >= 2017, ]


#put "NULL" in all blank cells in the df
#df <- df %>% mutate(across(where(is.character), ~replace(., is.na(.), "NULL")))
print(head(df))

df$CORRECTED_GENOTYPE <- NA #create a corrected genotype column that will be added to later

# these next several lines of code are to delete spaces within a single gene. This is necessary because I'm 
#using spaces to split the genes and if a single gene contains a space it will be split into multiple genes


# Function to process gene names
process_gene_names <- function(text) {
  # Step 1: Join parentheses with preceding word
  text <- gsub("(\\w+)\\s+(\\([^)]+\\))", "\\1\\2", text)
  
  # Step 2: Handle cases like "GyrA (S83L), (D87N)"
  text <- gsub("(\\w+\\([^)]+\\)),\\s*(\\([^)]+\\))", "\\1 \\1\\2", text)
  
  # Step 3: Split combined genes into separate genes
  text <- gsub("(\\w+\\([^)]+\\))\\([^)]+\\)", "\\1 \\1\\2", text)
  
  # Step 4: Remove any remaining commas
  text <- gsub(",", "", text)
  
  return(text)
}

# Apply the function to relevant columns
gene_columns <- colnames(df)[which(colnames(df) == "Aminoglycoside_Resist_Genes"):which(colnames(df) == "Other_Classes_Resist_Genes")]

for (col in gene_columns) {
  df[, col] <- sapply(df[, col], process_gene_names)
}

# Apply to GENOTYPE column as well
df[, "GENOTYPE"] <- sapply(df[, "GENOTYPE"], process_gene_names)



#sometimes genes that are classified elsewhere (i.e., for other isolates in the dataset) are not classified for a
#particular isolate (i.e., the gene only appears in the GENOTYPE column)
#in this case, I do not want to add that gene to the other_resistant_genes column
#these next few lines of code get a vector of all the classified genes so that I can check to see if a gene
#should be added to the other_resistant_genes column or not
not_other_genes_list <- list()
for (row in 1:nrow(df)) {
not_other_genes <- str_split(df[row, which(colnames(df) == "Aminoglycoside_Resist_Genes"):which(colnames(df) == "Phenicol_Resist_Genes")], " ")
not_other_genes_list <- list.append(not_other_genes_list, not_other_genes)
}
not_other_genes_vec <- unique(unlist(not_other_genes_list))
    

for (row in 1:nrow(df)) {

    #get all the genes from each of the gene class columns
    genes_vec <- unlist(str_split(df[row, which(colnames(df) == "Aminoglycoside_Resist_Genes"):which(colnames(df) == "Other_Classes_Resist_Genes")], " "))
    
    #get all the genes from the GENOTYPE column
    genotype_vec <- unlist(str_split(df[row, "GENOTYPE"], " "))

    #find which genes need to go into the other resistant genes column
    diff <- paste(unique(genotype_vec[! genotype_vec %in% genes_vec]), collapse = " ")
    if (! diff %in% not_other_genes_vec) {
    df[row, "Other_Classes_Resist_Genes"] <- diff
    }

    #get all the genes that need to go into the corrected genotype column
    unique_vec <- unique(c(genes_vec, genotype_vec))
    if (length(unique_vec) > 1) {
        unique_vec<- unique_vec[!(unique_vec == "NULL")] #delete null if there is at least one gene present
        }

    #add genes to the corrected genotype column
    unique_vec <- paste(unique_vec, collapse = " ")
    df[row, "CORRECTED_GENOTYPE"] <- unique_vec
}
if (FALSE) {
# Remove "NA" from the corrected genotype column unless that is the only value present
df$CORRECTED_GENOTYPE <- sapply(df$CORRECTED_GENOTYPE, function(x) {
  genes <- unlist(strsplit(x, " "))
  genes <- genes[genes != ""]  # Remove empty strings
  if (length(genes) == 1 && genes == "NA") {
    return("NA")
  } else {
    genes <- genes[genes != "NA"]
    return(paste(genes, collapse = " "))
  }
})
}



print(head(df$CORRECTED_GENOTYPE))
write.csv(df, file = str_glue("{data_source}/{data_source}_corrected_genotype.csv"), row.names = FALSE )
return(df)

}


#retail_meats_df <- read.csv("CVM-2020-NARMS-RetailMeatsData.csv")
#fix_genotype(retail_meats_df, "retail_meats")






MIC_to_interval <- function(df, data_source) {

    #only get the data we want to analyze
    df <- df[df$GROWTH == "YES", ]
    df <- df[df$GENUS == "EC", ]
    df <- df[df$SPECIES == "coli", ]
    df <- df[df$HOST_SPECIES == "Cattle", ]
    df <- df[df$Year != 2020, ]
    #only get 2017 onward
    df <- df[df$Year >= 2017, ]

    print(colnames(df))
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


    write.csv(merged_df, file = str_glue("{data_source}/{data_source}_phenotype_data.csv"), row.names = FALSE)

    return(merged_df)

}




prevalence_descriptives <- function(df, data_source, resistance_indicator, class_level = FALSE) {

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

    print("here1")
    prevalence_df <- bind_rows(df_list)
    print("here2")

    prevalence_df <- prevalence_df * 100
    prevalence_df <- cbind(N = unlist(sample_size_list), prevalence_df)
    print("here3")

    prevalence_df <- cbind(Year = unlist(year_list), prevalence_df)
    print("here4")
    if (class_level == FALSE) {
        write.csv(prevalence_df, str_glue("{data_source}/prevalence_descriptives/{data_source}_{resistance_indicator}_prevalence_descriptives.csv"), row.names = FALSE)
    } else if (class_level == TRUE) {
        write.csv(prevalence_df, str_glue("{data_source}/prevalence_descriptives/{data_source}_{resistance_indicator}_class_level_prevalence_descriptives.csv"), row.names = FALSE)
    }

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

    write.csv(prevalence_df, str_glue("{data_source}/prevalence_descriptives/{data_source}_{resistance_indicator}_class_level_prevalence_descriptives.csv"), row.names = FALSE)
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
            mutate(!!col_name := pmax(!!!rlang::syms(cols), na.rm = TRUE)) %>%
            select(-all_of(cols), !!col_name)
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
match_to_refgenes <- function(narms_df, refgenes_df) {
#get rid of rows where the Genes column is duplicated
narms_df <- narms_df[!duplicated(narms_df$Genes), ]
  
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
#narms_df <- read.csv("gene_class_mappings.csv")
#refgenes_df <- read.csv("refgenes.csv")

#match_to_refgenes(narms_df, refgenes_df)




convert_gene_to_gene_family <- function(df, genotype_class_mappings, data_source) {
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

  # Function to create a consistent representation for matching
  create_matching_name <- function(name) {
    gsub("[^A-Za-z0-9]", "", tolower(name))
  }

  # Create matching names for both the dataframe columns and the mapping file
  df_matching_names <- create_matching_name(names(df))
  mapping_matching_names <- create_matching_name(genotype_class_mappings$Genes)

  # Create a named vector for renaming, using original names but matching based on the consistent representation
  genotype_class_names <- setNames(
    make.unique(genotype_class_mappings$Gene_family),
    names(df)[match(mapping_matching_names, df_matching_names)]
  )

  print("Gene mappings:")
  print(genotype_class_names)

  genotype_duplicated_classes <- names(which(table(genotype_class_mappings$Gene_family) > 1))

  # Function to merge duplicated columns keeping the largest value
  merge_columns <- function(df, cols) {
    col_name <- unique(gsub("\\.\\d+$", "", cols[1]))
    df[[col_name]] <- apply(df[, cols, drop = FALSE], 1, max, na.rm = TRUE)
    df <- df[, !(names(df) %in% cols)]
    df[[col_name]] <- df[[col_name]]  # Add the new column
    return(df)
  }


  #merge any columns in df that are currently duplicated
  for (col in colnames(df)) {
    cols_to_merge <- names(df)[grepl(paste0("^", col, "(\\.|$)"), names(df))]
    
    if (length(cols_to_merge) > 1) {
      df <- merge_columns(df, cols_to_merge)
    }
  }

  # Rename and merge columns
  for (old_name in names(genotype_class_names)) {
    new_name <- genotype_class_names[old_name]
    if (old_name %in% names(df)) {
      names(df)[names(df) == old_name] <- new_name
    }
  }
  
  for (class in genotype_duplicated_classes) {
    cols_to_merge <- names(df)[grepl(paste0("^", class, "(\\.|$)"), names(df))]
    if (length(cols_to_merge) > 1) {
      df <- merge_columns(df, cols_to_merge)
    }
  }

  # Remove trailing dot and number from column names in genotype_df
  colnames(df) <- gsub("\\.\\d+$", "", colnames(df))

  # Check for unmapped genes
  unmapped_genes <- setdiff(names(df), c("ID", "Year", unique(genotype_class_mappings$Gene_family)))
  if (length(unmapped_genes) > 0) {
    warning("The following genes were not mapped to any class: ", paste(unmapped_genes, collapse = ", "))
  }

  # Debugging: Print the final column names
  print("Final column names:")
  print(colnames(df))
  write.csv(df, file = str_glue("{data_source}/{data_source}_wide_family_level_genotype.csv"), row.names = FALSE)
  
  return(df)
}



convert_gene_to_class <- function(df, genotype_class_mappings, data_source) {
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

 #genotype_class_mappings <- read.csv("gene_class_mappings.csv")
        genotype_class_mappings$Genes <- gsub(pattern = "[^A-Za-z0-9]", replacement = ".", genotype_class_mappings$Genes)
        #get rid of duplicate duplicates in Genes column
        genotype_class_mappings <- genotype_class_mappings[!duplicated(genotype_class_mappings$Genes), ]

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

      #merge any columns in df that are currently duplicated
      for (col in colnames(df)) {
                    cols_to_merge <- names(df)[grepl(paste0("^", col, "(\\.|$)"), names(df))]
                    
                    if (length(cols_to_merge) > 1) {
                        df <- merge_columns(df, cols_to_merge)
                    }
                }


        df <- df %>%
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
        colnames(df) <- gsub("\\.\\d+$", "", colnames(df))

        # Check for unmapped genes
        unmapped_genes <- setdiff(names(df), c("ID", "Year", unique(genotype_class_mappings$corrected_class)))
        if (length(unmapped_genes) > 0) {
            warning("The following genes were not mapped to any class: ", paste(unmapped_genes, collapse = ", "))
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
  write.csv(df, file = str_glue("{data_source}/{data_source}_wide_class_level_genotype.csv"), row.names = FALSE)
  
  return(df)
}



#cecal_df <- read.csv("cecal/cecal_wide_corrected_genotype.csv")
#retail_meats_df <- read.csv("Retail_Meats/Retail_Meats_wide_corrected_genotype.csv")
#genotype_class_mappings <- read.csv("gene_class_mappings.csv")
#convert_gene_to_class(cecal_df, genotype_class_mappings, data_source = "cecal")
#convert_gene_to_class(retail_meats_df, genotype_class_mappings, data_source = "Retail_Meats")






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
  write.csv(df, file = str_glue("{data_source}/{data_source}_wide_class_level_phenotype.csv"), row.names = FALSE)
  
  return(df)
}


# cecal_phenotype_df <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_resStatus_phenotype.csv")
# retail_meats_phenotype_df <- read.csv("Retail_Meats_2017Onward/Retail_Meats_2017Onward_wide_resStatus_phenotype.csv")
# phenotype_class_mappings <- read_excel("EcoliBreakPoints.xlsx")
# convert_phenotype_to_class(cecal_phenotype_df, phenotype_class_mappings, data_source = "cecal_2017Onward")
# convert_phenotype_to_class(retail_meats_phenotype_df, phenotype_class_mappings, data_source = "Retail_Meats_2017Onward")
