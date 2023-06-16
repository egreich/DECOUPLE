
# README for DECOUPLE Repo

Repository for scripts related to the FLUXNET Secondment project: Evaluating SIF-T and SIF-GPP Decoupling Across Timescales and Ecosystems. Data files accompanying this repo can be found here (tbd).

### Folder descriptions

* data_raw
  - messy data compiled from various PIs
* data_formatted
  - formatted/ organized data. Data from data_raw after running 00_format_data.R
* data_clean
  - Csv file data filtered for noise and quality issues. This is the data from data_formatted after running 01_prep_data.R.
* scripts
* shell_scripts
  - bash scripts for running models on an HPC (slurm).
* END files and Slurm_jobs_3chains.csv
  - files describing parallel processes run using an HPC. These files are called in the shell scripts
* models
  - JAGS scripts for each model, saved as .R files
* output_ folders
  - model outputs
* plots
