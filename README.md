# Summary
This repository contains code for the analyses conducted in "Clinical features, genetics, and pathology in a large series of movement disorder cases: a retrospective multi-ancestry brain bank cohort study".
The aim of the study is to assess clinico-pathological correlation in individuals carrying disease-associated and risk genetic variants, the frequency of clinical misdiagnosis and the impact of co-pathology.

## Data Statement 
- All MDGAP/GP2 data are hosted in collaboration with the Accelerating Medicines Partnership in Parkinson's Disease and are available via application on the website.
- The clinical and genetic data is available via the GP2 website (https://gp2.org). Genotyping imputation, quality control, ancestry prediction, and processing were performed using GenoTools (v1.3.5), publicly available on GitHub.
- Pathological staging and co-pathology data will be included in Release 11 of GP2, scheduled for December 2025.

## Helpful Links
- https://gp2.org/
- https://amp-pd.org/

# Repository Orientation
- The analyses/ directory includes all analyses discussed in the manuscript
<pre> THIS_REPO/ 
  ├── analyses/ 
  |     ├── 00_Data_Clean.R
  |     └── 01_Sensitivity_Specificity.R
  |     └── 02_Copathologies.R
  |     └── 03_Genes_Extraction.R
  |     └── 04_Path_Stages.R
  |     └── 05_Ancestry.R
  └── README.md 

  # Software
  | Software |        Version(s)    |     Resource URL     |     RRID    |      Notes    |
|-----------|----------------------|---------------------|--------------|----------------|
|R Project for Statistical Computing| 4.5.1  | (http://www.r-project.org/)| RRID:SCR_001905 | tidyverse; dplyr; tidyr; ggplot; data.table; used for general data wrangling/plotting/analyses |
|PLINK| 2.0  |  (https://www.cog-genomics.org/plink/2.0/) | RRID:SCR_001757 | used for genetic analyses |
|samtools (bcftools)| 1.9 |  (https://samtools.github.io/bcftools/) | RRID:SCR_002105	 | VCF maniplulation |

