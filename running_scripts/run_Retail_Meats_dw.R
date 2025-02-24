library(here)

source(here("function_scripts", "data_wrangling.R"))
source(here("function_scripts", "data_wrangling_Retail.R"))


#RUN THE DATA WRANGLING FOR THE PHENOTYPE
Retail_Meats_df <- read.csv("Retail_Meats_Data/CVM-NARMS-Retail-Meats.csv")

MIC_df <- Retail_Meats_MIC_to_interval(df = Retail_Meats_df, data_source = "Retail_Meats")
breakpoint_df <- add_break_points(df = MIC_df, data_source = "Retail_Meats")
phenotype_wide_df <- add_resistance_status(df = breakpoint_df, data_source = "Retail_Meats")

phenotype_class_mappings <- read_excel("EcoliBreakPoints.xlsx")
convert_phenotype_to_class(df = phenotype_wide_df, phenotype_class_mappings, data_source = "Retail_Meats")







#RUN THE DATA WRANGLING FOR THE GENOTYPE
gene_df <- read.csv("Gene_Data/amrfinder_results-AMR-CORE.point-muts-included.csv")
ID_Year_lookup_df <- read.csv("Retail_Meats_Data/CVM-NARMS-Retail-Meats.csv")
Retail_Meats_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "Retail_Meats", output_level = "class")
Retail_Meats_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "Retail_Meats", output_level = "subclass")

wide_gene_df <- Retail_Meats_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "Retail_Meats", output_level = "gene")

#get gene family
genotype_class_mappings <- read.csv("gene_class_mappings.csv")

convert_gene_to_gene_family(df = wide_gene_df, genotype_class_mappings, data_source = "Retail_Meats")






#####
#DO 2017 ONWARD TOO
#####


#RUN THE DATA WRANGLING FOR THE PHENOTYPE
MIC_df2017 <- Retail_Meats_MIC_to_interval(df = Retail_Meats_df, data_source = "Retail_Meats_2017Onward", limit_date_range = TRUE)
breakpoint_df2017 <- add_break_points(df = MIC_df2017, data_source = "Retail_Meats")
phenotype_wide_df2017 <- add_resistance_status(df = breakpoint_df2017, data_source = "Retail_Meats_2017Onward")

convert_phenotype_to_class(df = phenotype_wide_df2017, phenotype_class_mappings, data_source = "Retail_Meats_2017Onward")







# RUN THE DATA WRANGLING FOR THE GENOTYPE
Retail_Meats_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "Retail_Meats_2017Onward", output_level = "class", limit_date_range = TRUE)
Retail_Meats_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "Retail_Meats_2017Onward", output_level = "subclass")

wide_gene_df_2017 <- Retail_Meats_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "Retail_Meats_2017Onward", output_level = "gene", limit_date_range = TRUE)

# get gene family
convert_gene_to_gene_family(df = wide_gene_df_2017, genotype_class_mappings, data_source = "Retail_Meats_2017Onward")

