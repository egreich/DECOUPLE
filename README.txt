
Repository for scripts related to the FLUXNET Secondment project: Evaluating SIF-T and SIF-GPP Decoupling Across Timescales and Ecosystems. Data files accompanying this repo can be found here (tbd).

## Folder descriptions
* data_raw
- messy data compiled from various PIs
* data_formatted
- formatted/ organized data. Data after running 00_format_data.R
* data_clean
- Csv file data filtered for noise and quality issues. This is the formatted data after running 01_filter_data.R. These files are used as inputs to run 02_prep_data. The output data from 02_prep_data are saved in this folder as RData files.
* scripts
* shell_scripts
- bash scripts for running models on an HPC (slurm).
* models
* plots

