library(here)

source(here("function_scripts", "data_wrangling.R"))

#RUN THE DATA WRANGLING FOR THE PHENOTYPE
phenotype_wide_df <- read.csv("cecal/cecal_wide_resStatus_phenotype.csv")
phenotype_class_mappings <- read_excel("EcoliBreakPoints.xlsx")
convert_phenotype_to_class(df = phenotype_wide_df, phenotype_class_mappings, data_source = "cecal")

#RUN THE DATA WRANGLING FOR THE GENOTYPE
wide_gene_df <- read.csv("cecal/cecal_wide_corrected_genotype.csv")
genotype_class_mappings <- read.csv("gene_class_mappings.csv")

#get gene family
convert_gene_to_gene_family(df = wide_gene_df, genotype_class_mappings, data_source = "cecal")

#get gene classes
convert_gene_to_class(df = wide_gene_df, genotype_class_mappings, data_source = "cecal")




#####
#DO 2017 ONWARD TOO
#####


#RUN THE DATA WRANGLING FOR THE PHENOTYPE
phenotype_wide_df2017 <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_resStatus_phenotype.csv")
phenotype_class_mappings <- read_excel("EcoliBreakPoints.xlsx")
convert_phenotype_to_class(df = phenotype_wide_df2017, phenotype_class_mappings, data_source = "cecal_2017Onward")

# #RUN THE DATA WRANGLING FOR THE GENOTYPE
wide_gene_df2017 <- read.csv("cecal_2017Onward/cecal_2017Onward_wide_corrected_genotype.csv")
genotype_class_mappings <- read.csv("gene_class_mappings.csv")

# #get gene family
convert_gene_to_gene_family(df = wide_gene_df2017, genotype_class_mappings, data_source = "cecal_2017Onward")

#get gene classes
convert_gene_to_class(df = wide_gene_df2017, genotype_class_mappings, data_source = "cecal_2017Onward")
