library(here)

source(here("function_scripts", "data_wrangling_NAHLN.R"))
source(here("function_scripts", "data_wrangling.R"))

#RUN THE DATA WRANGLING FOR THE PHENOTYPE
phenotype_wide_df <- wrangle_NAHLN_phenotype()
phenotype_class_mappings <- read_excel("EcoliBreakPointsWithEcoffs.xlsx")
convert_phenotype_to_class(df = phenotype_wide_df, phenotype_class_mappings, data_source = "NAHLN")

#RUN THE DATA WRANGLING FOR THE GENOTYPE
wide_gene_df <- wrangle_NAHLN_genotype()
genotype_class_mappings <- read.csv("wNAHLN_gene_class_mappings.csv")

#get gene family
convert_gene_to_gene_family(df = wide_gene_df, genotype_class_mappings, data_source = "NAHLN")

#get gene classes
convert_gene_to_class(df = wide_gene_df, genotype_class_mappings, data_source = "NAHLN")