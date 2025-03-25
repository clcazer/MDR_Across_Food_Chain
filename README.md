# Analyzing multidrug resistance patterns across the food supply chain using association rule mining
Association mining on MDR patterns in cattle-associated E. coli.

__Abstract__
We used the machine learning method association rule mining to analyze multidrug resistance (MDR) among   cattle-associated Escherichia coli along the food supply chain in the USA. All datasets were stratified by year, source, and resistance indicator (genotypic/phenotypic). Pruned rulesets were compared by calculating the proportion of rules from a comparison ruleset that are captured in a reference ruleset. Rulesets were compared across years within each source and indicator type to quantify how MDR patterns change over time. At the class level, on average nearly 50% or more of the MDR patterns remain the same year over year for genotypic and phenotypic indicators. Rulesets were compared between data sources to quantify how MDR patterns change across the food supply chain. These comparisons suggest that there is a greater diversity of MDR patterns present at slaughterhouse settings than at retail settings; and further, that there is a greater diversity of MDR patterns amongst sick cattle on farm settings than at either slaughterhouse or retail settings. Genetic evidence supports this being attributable to a greater genetic diversity associated with pathogenic bacteria vs commensals. Rulesets were compared between indicators to quantify the degree of correspondence between phenotypic and genotypic data. Genotypic rulesets were better able to capture phenotypic rulesets than the reverse. Adding another aminoglycoside (streptomycin) to the phenotypic analysis, improved ruleset correspondence. This asymmetry may be driven by drug specific aminoglycoside resistance genes, suggesting that more drugs need to be assessed to have a fuller understanding of the variation in MDR patterns.

Joshua Glass[1#], Gayatri Anil[2], Jordan D. Zehr[3], Kristina M. Ceres[4], Laura B. Goodman[5], Casey L. Cazer[6]

1. #Joshua Glass, M.S., Department of Clinical Sciences, Cornell University College of Veterinary Medicine, Ithaca, NY 14850, United States, jg2527@cornell.edu
2. Gayatri Anil, Department of Clinical Sciences, Cornell University College of Veterinary Medicine, Ithaca, NY 14850, United States, ga325@cornell.edu
3. Jordan D. Zehr, Ph.D., Department of Public and Ecosystem Health, College of Veterinary Medicine, Cornell University, Ithaca, NY 14850, United States, jdz34@cornell.edu
4. Kristina M. Ceres, Ph.D., Department of Public and Ecosystem Health, College of Veterinary Medicine, Cornell University, Ithaca, NY 14850, United States, kc649@cornell.edu
5. Laura B. Goodman, Ph.D., Department of Public and Ecosystem Health, College of Veterinary Medicine, Cornell University, Ithaca, NY 14850, United States, laura.goodman@cornell.edu
6. Casey L. Cazer, DVM, Ph.D., Department of Clinical Sciences, Cornell University College of Veterinary Medicine, Ithaca, NY 14850, United States, casey.cazer@cornell.edu

\# Corresponding Author


## Table of Contents
- [Data Description and Availability](#data-description-and-availability)
- [Environment Setup](#environment-setup)
- [Repository Organization](#repository-organization)
- [Steps to Reproduce Analysis](#steps-to-reproduce-analysis)
- [Licencse](#license)


## Data Description and Availability
Data from two antimicrobial surveilance organizations were used. One organization was The National Antimicrobial Resistance Monitoring System [(NARMS)](https://www.fda.gov/media/79976/download?attachment). The NARMS program is a national antimicrobial resistance surveillance program collecting samples from slaughterhouses and retail settings, focused on foodborne pathogens. This data is *publicly available* and can be found in the ```Slaughterhouse_Data``` folder and the ```Retail_Meats_Data``` folder, as well as here:

>Food and Drug Administration (FDA), 2024. NARMS Now [WWW Document]. Rockville, MD: U.S. Department of Health and Human Services. Available from URL: https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/narms-now-integrated-data. Accessed 12/03/2024.

The second organization was National Animal Health Laboratory Network [(NAHLN)](https://www.aphis.usda.gov/sites/default/files/2020_APHIS_AMR_PilotProject_EOY_Report.pdf). The NAHLN antimicrobial resistance program collects data from participating veterinary diagnostic laboratories on bacteria isolated from animal samples submitted as clinical cases. __This data was used under a *private data sharing agreement* and cannot be shared publicly.__

The data consists of bacterial isolates tested against several antimicrobials. The analysis only examines *Escherichia coli* isolates associated with cattle, the filtering for this is accomplished in a function described later. Genotypic data corresponding to a subset of isolates was also analyzed, the genotypic data can be found in the ```Gene_Data``` folder.

The data is divided into three datasets:
1. Retail Meats: this is data from NARMS, consisting of samples taken from retail meats. This dataset is referred to as "Retail_Meats" throughout the code.
2. Slaughterhouse_Data: this is data from NARMS, consisting of cecal samples taken from animal at slaughter. This dataset is referred to as "cecal" thoughout the code.
3. NAHLN Data: this is data from NAHLN, consisting of samples taken from clinically ill animals. This is referred to as "NAHLN" throughout the code, but is referred to as the "sick cattle" dataset in the paper.

## Environment Setup
This project uses renv for package management. To set up the environment:

1. Install R and RStudio
2. Install renv:
```install.packages("renv")```
3. Clone this repository
4. Open the project in RStudio
5. Run:
```renv::restore()```

This will automatically install all required packages with the exact versions used in the analysis.

## Repository Organization


```
├── Raw_Data/ # Raw phenotypic and genotypic data
│ ├── Gene_Data/ # AMR gene data files
│ ├── NAHLN_Data/ # Clinical isolate data (not included due to data privacy)
│ ├── Retail_Meats_Data/ # Retail meat sample data
│ ├──  Slaughterhouse_Data/ # Cecal sample data
│ └── Virulence_Plasmid_data/ # data used for the virulence gene and plasmid analyses
|
├── Dataset_Specific_Outputs/ # Analysis outputs by dataset
│ ├── cecal/ # Slaughterhouse data analysis
│ │ ├── figures/ # Visualization outputs
│ │ ├── prevalence_descriptives/ # Statistical summaries
│ │ ├── ruleData/ # Mining results
│ │ └── tabulation_comp/ # Comparative analyses
│ ├── Retail_Meats/ # Similar structure for retail data
│ ├── NAHLN/ # Similar structure for clinical data
│ ├── cecal_2017Onward/ # Year>=2017 slaughterhouse analysis
│ └── Retail_Meats_2017Onward/ # Year>=2017 retail analysis
│
├── combined_outputs/ # Cross-dataset analyses (i.e., the outputs used in the manuscript)
│ ├── combined_csv_outputs/ # Consolidated csv results (e.g., descriptive stats tables; consolidated manually)
│ └── combined_figures_outputs/ # Consolidated figure outputs (where the two figures for the mauscripts are saved)
│
├── Scripts/ # Analysis code
│ ├── function_scripts/ # Core functions and utilities
│ └── running_scripts/ # Analysis execution scripts
│
└── Virulence_Plasmid_outputs/ # Virulence and plasmid analyses
 ├── virulence_descriptives/ # Virulence gene statistics
 └── plasmid_rule_comp/ # Plasmid pattern analysis
```





#### Data Directories
- ```Raw_Data``` contains all shareable phenotypic and genotypic data.

#### Dataset Specific Output Directories
- ```cecal```
- ```Retail_Meats```
- ```NAHLN```

_Cecal and retail meats datasets included the drug streptomycin; for which breakpoints changed in 2014. The orgininal analysis excludes streptomycin due incompatabilities resulting analyzing data across the time of the breakpoint change. However, genotypic data for both of these datasets begins in 2017 (after the breakpoints changed) so streptomycin could be added back in without problems in order to better examine phenotypic and genotypic correspondance. For this purpose the following two directoies exist:_

- ```cecal_2017Onward```
- ```Retail_Meats_2017Onward```

#### Combined Output Directories
These are where all the combined outputs that are used in the manuscript can be found:
- ```combined_csv_outputs```
- ```combined_figures_outputs```

#### Script Directories

- ```function_scripts``` contains files several files. The files contain several functions. Each file is named based on the type of function it contains, such as ```data_wrangling.R``` or ```rule_mining_and_selection.R```
- ```running_scripts``` contains scripts for running the functions that are defined in the files from the ```function_scripts``` directory. There are several different running scripts that accomplish different portions of the analysis- they are named for the task that they accomplish. For example, ```run_Retail_Meats_dw.R``` runs all the data wrangling functions necessary for the retail meats dataset, while ```run_Retail_Meats.R``` runs the actual analysis for the pre-processed retail meats dataset.

### Virulence and Plasmid Outputs
This directory contains outputs for an analysis comparing the counts of virulence genes across each dataset, as well as outputs for analysis comparing gene rule combinations to gene combinations found on plasmids.

## Steps to Reproduce Analysis

1. __Run Data Wrangling Code__
    - execute ```Scripts\running_scripts\run_Cecal_dw.R```
    - execute ```Scripts\running_scripts\run_Retail_Meats_dw.R```
    - execute ```Scripts\running_scripts\run_NAHLN_dw.R```

2. __Run Dataset Specific Analysis Code:__ this will produce all of the dataset specific outputs including the descriptive stats.
    - execute ```Scripts\running_scripts\run_Cecal.R```
    - execute ```Scripts\running_scripts\run_Cecal_2017On.R```
    - execute ```Scripts\running_scripts\run_Retail_Meats.R```
    - execute ```Scripts\running_scripts\run_Retail_Meats_2017On.R```
    - execute ```Scripts\running_scripts\run_NAHLN.R```

3. __Run Combined Outputs:__ this will produce all of the combined outputs that are used in the manuscript.
    - execute ```Scripts\running_scripts\run_Percent_Capture.R```
    - execute ```Scripts\function_scripts\make_graph_panel.R```
    - execute ```Scripts\function_scripts\tabulation_comp.R```
    - retreive all descriptive stats from the dataset specific directories and aggregate them manually

4. __Run Virulence Gene and Plasmid Analysis__
    - execute ```Scripts\function_scripts\virulence_analysis.R```
    - execute ```Scripts\function_scripts\plasmid_analysis.R```



## License
BSD 3-Clause License

Copyright (c) 2025, Joshua Glass

All rights reserved.


Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


