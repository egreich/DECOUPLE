
# README for DECOUPLE Repo

Repository for scripts related to the FLUXNET Secondment project: Decoupling of transpiration, gross primary productivity, and solar-induced fluorescence at the ecosystem scale.

### Folder descriptions

* scripts
  - script workflow labeled as 00a-02b, from data cleaning to running the models
* shell_scripts
  - bash scripts for running models on an HPC (slurm).
* END files and Slurm_jobs_3chains.csv
  - files describing parallel processes run using an HPC. These files are called in the shell scripts
* models
  - JAGS scripts for each model, saved as .R files
* output_ folders
  - model outputs
