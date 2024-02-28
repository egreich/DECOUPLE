#!/bin/bash
#SBATCH --job-name=mod2_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/DECOUPLE/log/Mod2_%A_%a.log
#SBATCH --cpus-per-task=3
#SBATCH --time=60:00:00
#SBATCH --mem=40000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=13-22
##SBATCH --array=1-24

### %A is monsoon job number %a is interior array index

module load R/4.1.2 # load a specific R version

chmod +x shell_scripts/mod2_job.sh # for permissions
chmod +x scripts/02b_run_model2.R # for permissions

##site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteEND)
varname=$(sed -n "$SLURM_ARRAY_TASK_ID"p mod2varnameEND)
sitename=$(sed -n "$SLURM_ARRAY_TASK_ID"p mod2sitenameEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p mod2seedEND)
scale=$(sed -n "$SLURM_ARRAY_TASK_ID"p mod2scaleEND)

# Run the analysis
srun ./shell_scripts/mod2_job.sh $varname $sitename $seed $scale
