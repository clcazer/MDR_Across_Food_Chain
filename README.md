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
KJAHSFKJSBFKFJED

### Description of Directories
### Description of Files
### Description of Functions


## Steps to Reproduce Analysis