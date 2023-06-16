#!/bin/bash
#SBATCH --job-name=mod1_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/DECOUPLE/log/AllMod1_%A_%a.log
#SBATCH --cpus-per-task=3
#SBATCH --time=6:00:00
#SBATCH --mem=80000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=1-35

### %A is monsoon job number %a is interior array index

module load R/4.1.2 # load a specific R version

chmod +x shell_scripts/mod1_job.sh # for permissions
chmod +x scripts/02_run_model.R # for permissions

##site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteEND)
varname=$(sed -n "$SLURM_ARRAY_TASK_ID"p varnameEND)
sitename=$(sed -n "$SLURM_ARRAY_TASK_ID"p sitenameEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedEND)

# Run the analysis
srun ./shell_scripts/mod1_job.sh $varname $sitename $seed
