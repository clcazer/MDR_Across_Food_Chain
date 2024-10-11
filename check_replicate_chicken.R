library(stringr)
library(rlist)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(readxl)
library(cowplot)
library(arules)
library(dplyr)
source("data_wrangling.R")
source("rule_mining_and_selection.R")
source("rule_analysis.R")



MIC_to_interval <- function(df) {

    #only get the data we want to analyze
    #df <- df[df$GROWTH == "YES", ]
    df <- df[df$GENUS == "EC", ]
    df <- df[df$SPECIES == "coli", ]
    df <- df[df$HOST_SPECIES == "Chickens", ]
    df <- df[df$Year >= 2004, ]
    df <- df[df$Year <= 2012, ]

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
    "NAL.Sign",
    "SMX.Sign",
    "TET.Sign",
    

    "TIO.Sign",
    "AMI.Sign",
    "STR.Sign",
    "KAN.Sign"
    )]



   # write.csv(df, file = "df_test.csv", row.names = FALSE)


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
    "NAL",
    "TET",

    "TIO",
    "AMI",
    "Str",
    "KAN"
    )]

    #pivot the sign_df to long format delete the .Sign from the Drugcode column so we are just left with drug code that will match that on the MIC_df
    sign_df <- pivot_longer(sign_df, cols = ends_with("Sign"), names_to = "Drugcode", values_to = c("MIC_SIGN"), values_drop_na = FALSE)# pivot to long format
    sign_df$Drugcode  <- substr(sign_df$Drugcode , 0, 3)
    sign_df$Drugcode <- toupper(sign_df$Drugcode) # make everything uppercase


    #pivot the MIC_df to long format
    MIC_df <- pivot_longer(MIC_df, cols = -c(ID, Year, SOURCE), names_to = "Drugcode", values_to = c("MIC"), values_drop_na = TRUE)# pivot to long format
    MIC_df$Drugcode <- toupper(MIC_df$Drugcode) # make everything uppercase

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

    print(length(unique(merged_df$ID)))
    write.csv(merged_df, file = "retail_chicken_ecoli_phenotype_data.csv", row.names = FALSE)

}




#df <- read.csv("CVM-2020-NARMS-RetailMeatsData.csv")
#MIC_to_interval(df = df)


#df <- read.csv("Retail_Meats_Chicken_Replicate/retail_chicken_ecoli_phenotype_data.csv")
#new_df <- add_break_points(df = df, data_source = "Retail_Meats_Chicken_Replicate")

#add_resistance_status(df = new_df, data_source = "Retail_Meats_Chicken_Replicate")


df <- read.csv("Retail_Meats_Chicken_Replicate/Retail_Meats_Chicken_Replicate_wide_resStatus_phenotype_NOTdb.csv")


#plot_all_vs_best_rules(df = df, resistance_indicator = "phenotype", 
                        #target = "rules", cut_off = c(0.75, 2, 0.5, 0.01),
                        #measures_used = c("confidence", "lift", "phi", "support"), data_source = "Retail_Meats_Chicken_Replicate")



#rule_overlap(df = df, resistance_indicator = "phenotype", 
                       # target = "rules", cut_off = c(0.75, 2, 0.5, 0.01),
                        #measures_used = c("confidence", "lift", "phi", "support"), data_source = "Retail_Meats_Chicken_Replicate", rules_selected = "best")


#cumulative_rule_stability(df = df, resistance_indicator = "phenotype", 
                        #target = "rules", cut_off = c(0.75, 2, 0.5, 0.01),
                        #measures_used = c("confidence", "lift", "phi", "support"), data_source = "Retail_Meats_Chicken_Replicate", rules_selected = "best")




#prevalence_descriptives(df = df, data_source = "Retail_Meats_Chicken_Replicate")