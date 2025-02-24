library(here)

source(here("function_scripts", "data_wrangling_NAHLN.R"))
source(here("function_scripts", "data_wrangling.R"))

#RUN THE DATA WRANGLING FOR THE PHENOTYPE
df <- read_excel("NAHLN_Data/NAHLN_CattleEcoli_PPY1-PPY5_10-12-2023.xlsx")
break_points_df <- read_excel("EcoliBreakPointsWithEcoffs.xlsx")
phenotype_wide_df <- wrangle_NAHLN_phenotype(df, break_points_df)


phenotype_class_mappings <- read_excel("EcoliBreakPointsWithEcoffs.xlsx")
convert_phenotype_to_class(df = phenotype_wide_df, phenotype_class_mappings, data_source = "NAHLN")





#RUN THE DATA WRANGLING FOR THE GENOTYPE
gene_df <- read.csv("Gene_Data/amrfinder_results-AMR-CORE.point-muts-included.csv")
ID_lookup_df <- read_excel("NAHLN_Data/Isolate+BioSample.xlsx")
Year_lookup_df <- read.csv("NAHLN_Data/NAHLN_CattleEcoli_PPY1-PPY5_10-12-2023.csv")
NAHLN_process_genotype_data(gene_df = gene_df, Year_lookup_df = Year_lookup_df, ID_lookup_df = ID_lookup_df, data_source = "NAHLN", output_level = "class")
NAHLN_process_genotype_data(gene_df = gene_df, Year_lookup_df = Year_lookup_df, ID_lookup_df = ID_lookup_df, data_source = "NAHLN", output_level = "subclass")



wide_gene_df <- NAHLN_process_genotype_data(gene_df = gene_df, Year_lookup_df = Year_lookup_df, ID_lookup_df = ID_lookup_df, data_source = "NAHLN", output_level = "gene")



#get gene family
genotype_class_mappings <- read.csv("gene_class_mappings.csv")

convert_gene_to_gene_family(df = wide_gene_df, genotype_class_mappings, data_source = "NAHLN")

