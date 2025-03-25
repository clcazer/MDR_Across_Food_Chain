library(here)

source(here("Scripts", "function_scripts", "data_wrangling.R"))
source(here("Scripts", "function_scripts", "data_wrangling_cecal.R"))


#RUN THE DATA WRANGLING FOR THE PHENOTYPE
cecal_df <- read.csv("Raw_Data/Slaughterhouse_Data/CVM-NARMS-Cecal.csv")
MIC_df <- cecal_MIC_to_interval(df = cecal_df, data_source = "cecal")

breakpoint_df <- add_break_points(df = MIC_df, data_source = "cecal")
phenotype_wide_df <- add_resistance_status(df = breakpoint_df, data_source = "cecal")

phenotype_class_mappings <- read_excel("EcoliBreakPoints.xlsx")
convert_phenotype_to_class(df = phenotype_wide_df, phenotype_class_mappings, data_source = "cecal")









#RUN THE DATA WRANGLING FOR THE GENOTYPE
 gene_df <- read.csv("Raw_Data/Gene_Data/amrfinder_results-AMR-CORE.point-muts-included.csv")
 ID_Year_lookup_df <- read.csv("Raw_Data/Slaughterhouse_Data/old/CVM-2020-NARMS-Cecal-Data.csv")
 cecal_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "cecal", output_level = "class")
 cecal_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "cecal", output_level = "subclass")
 wide_gene_df <- cecal_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "cecal", output_level = "gene")


genotype_class_mappings <- read.csv("gene_class_mappings.csv")

#get gene family
convert_gene_to_gene_family(df = wide_gene_df, genotype_class_mappings, data_source = "cecal")






#####
#DO 2017 ONWARD TOO
#####


#RUN THE DATA WRANGLING FOR THE PHENOTYPE
MIC_df2017 <- cecal_MIC_to_interval(df = cecal_df, data_source = "cecal", limit_date_range = TRUE)
breakpoint_df2017 <- add_break_points(df = MIC_df2017, data_source = "cecal")
phenotype_wide_df2017 <- add_resistance_status(df = breakpoint_df2017, data_source = "cecal_2017Onward")

convert_phenotype_to_class(df = phenotype_wide_df2017, phenotype_class_mappings, data_source = "cecal_2017Onward")









# RUN THE DATA WRANGLING FOR THE GENOTYPE
 cecal_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "cecal_2017Onward", output_level = "class", limit_date_range = TRUE)
 cecal_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "cecal_2017Onward", output_level = "subclass", limit_date_range = TRUE)
 wide_gene_df_2017 <- cecal_process_genotype_data(gene_df, ID_Year_lookup_df, data_source = "cecal_2017Onward", output_level = "gene", limit_date_range = TRUE)

#get gene family
convert_gene_to_gene_family(df = wide_gene_df_2017, genotype_class_mappings, data_source = "cecal_2017Onward")


print("DONE RUNNING CECAL DATA WRANGLING")