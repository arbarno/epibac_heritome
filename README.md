# epibac_heritome
Scripts and important data for "Heritability of epigenome changes in response to probiotic inoculation in corals"

## Methylation analysis
This repository contains the scripts (or references to scripts) and input files used for analyzing the epigenetic response of *Acropora* sp. and its associate gametes to bacteria inoculation in the paper "Heritability of epigenome changes in response to probiotic inoculation in corals"

The pipeline to create the reference assembly used in this manuscript can be found here: https://github.com/arbarno/DNA-methylation

### The files in the main repository (and how they were obtained):
- `bash_methylation_analysis.sh`
  - Contains the bash pipeline used to extract the methylation info, resulting in the output files in the output/ folder
- `R_methylation_analysis.R`
  - Contains the R pipeline used to further analyze and visualize the results from the output files generated from the bash pipeline

### The folders in the main repository:
- `python_scripts`
  - The python scripts used at the end of the bash analysis
- `output`
  - Output files generated from the bash pipeline
- `ref_files`
  - Files related to the Acropora reference genome assembled for this study
