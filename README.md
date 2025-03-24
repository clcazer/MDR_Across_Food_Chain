# Analyzing multidrug resistance patterns across the food supply chain using association rule mining
Association mining on MDR patterns in cattle-associated E. coli.

## Table of Contents
- [Data Description and Availability](#data-description-and-availability)
- [Repository Organization](#repository-organization)
    - [Description of Directories](#description-of-directories)
    - [Description of Files](#description-of-files)
    - [Description of Functions](#description-of-functions)
- [Steps to Reproduce Analysis](#steps-to-reproduce-analysis)

## Data Description and Availability
Data from two antimicrobial surveilance organizations were used. One organization was The National Antimicrobial Resistance Monitoring System [(NARMS)](https://www.fda.gov/media/79976/download?attachment). The NARMS program is a national antimicrobial resistance surveillance program collecting samples from slaughterhouses and retail settings, focused on foodborne pathogens. This data is *publicly available* and can be found in the ```Slaughterhouse_Data``` folder and the ```Retail_Meats_Data``` folder, as well as here:

>Food and Drug Administration (FDA), 2024. NARMS Now [WWW Document]. Rockville, MD: U.S. Department of Health and Human Services. Available from URL: https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/narms-now-integrated-data. Accessed 12/03/2024.

The second organization was National Animal Health Laboratory Network [(NAHLN)](https://www.aphis.usda.gov/sites/default/files/2020_APHIS_AMR_PilotProject_EOY_Report.pdf). The NAHLN antimicrobial resistance program collects data from participating veterinary diagnostic laboratories on bacteria isolated from animal samples submitted as clinical cases. __This data was used under a *private data sharing agreement* and cannot be shared publicly.__

The data consists of bacterial isolates tested against several antimicrobials. The analysis only examines *Escherichia coli* isolates associated with cattle, the filtering for this is accomplished in a function described later. Genotypic data corresponding to a subset of isolates was also analyzed, the genotypic data can be found in the ```Gene_Data``` folder.

The data is divided into three datasets:
1. Retail Meats: this is data from NARMS, consisting of samples taken from retail meats. This dataset is referred to as "Retail_Meats" throughout the code.
2. Slaughterhouse_Data: this is data from NARMS, consisting of cecal samples taken from animal at slaughter. This dataset is referred to as "cecal" thoughout the code.
3. NAHLN Data: this is data from NAHLN, consisting of samples taken from clinically ill animals. This is referred to as "NAHLN" throughout the code, but is referred to as the "sick cattle" dataset in the paper.



## Repository Organization

### Description of Directories
This repository is organized into 11 main directories.

```
├── Raw_Data/ # Raw phenotypic and genotypic data
│ ├── Gene_Data/ # AMR gene data files
│ ├── NAHLN_Data/ # Clinical isolate data
│ ├── Retail_Meats_Data/ # Retail meat sample data
│ └── Slaughterhouse_Data/ # Cecal sample data
│
├── Dataset_Specific_Outputs/ # Analysis outputs by dataset
│ ├── cecal/ # Slaughterhouse data analysis
│ │ ├── figures/ # Visualization outputs
│ │ ├── prevalence_descriptives/ # Statistical summaries
│ │ ├── ruleData/ # Mining results
│ │ └── tabulation_comp/ # Comparative analyses
│ ├── Retail_Meats/ # Similar structure for retail data
│ └── NAHLN/ # Similar structure for clinical data
│
├── Time_Period_Specific/ # Analyses for specific time ranges
│ ├── cecal_2017Onward/ # Post-2017 slaughterhouse analysis
│ └── Retail_Meats_2017Onward/ # Post-2017 retail analysis
│
├── Combined_Outputs/ # Cross-dataset analyses
│ ├── combined_csv_outputs/ # Consolidated results
│ └── combined_percent_captured/ # Comparative metrics
│
├── Scripts/ # Analysis code
│ ├── function_scripts/ # Core functions and utilities
│ └── running_scripts/ # Analysis execution scripts
│
└── Virulence_Plasmid_data/ # Virulence and plasmid analyses
├── virulence_descriptives/ # Virulence gene statistics
└── plasmid_rule_comp/ # Plasmid pattern analysis
```





#### Data Directories
- ```Raw_Data``` contains all shareable phenotypic and genotypic data (other than virulence gene and plasmid info).
- ```Virulence_Plasmid_data``` contains data on vurulence genes and plasmid gene combinations. This directory also contains analysis outputs (within subdirectories ```Virulence_Plasmid_data\virulence_descriptives``` and ```Virulence_Plasmid_data\plasmid_rule_comp```)

#### Dataset Specific Output Directories
- ```cecal```
- ```Retail_Meats```
- ```NAHLN```

_Cecal and retail meats datasets included the drug streptomycin; for which breakpoints changed in 2014. The orgininal analysis excludes streptomycin due incompatabilities resulting analyzing data across the time of the breakpoint change. However, genotypic data for both of these datasets begins in 2017 (after the breakpoints changed) so streptomycin could be added back in without problems in order to better examine phenotypic and genotypic correspondance. For this purpose the following two directoies exist:_

- ```cecal_2017Onward```
- ```Retail_Meats_2017Onward```

#### Combined Output Directories

- ```combined_csv_outputs```
- ```combined_percent_captured```

#### Script Directories

- ```function_scripts``` contains files several files. The files contain several functions. Each file is named based on the type of function it contains, such as ```function_scripts\data_wrangling.R``` or ```function_scripts\rule_mining_and_selection.R```
- ```running_scripts``` contains scripts for running the functions that are defined in the files from the ```function_scripts``` directory. There are several different running scripts that accomplish different portions of the analysis- they are named for the task that they accomplish. For example, ```running_scripts\run_Retail_Meats_dw.R``` runs all the data wrangling functions necessary for the retail meats dataset, while ```running_scripts\run_Retail_Meats.R``` runs the actual analysis for the pre-processed retail meats dataset.


### Description of Files
### Description of Functions


## Steps to Reproduce Analysis