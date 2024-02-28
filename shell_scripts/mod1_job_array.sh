#!/bin/bash
#SBATCH --job-name=mod1_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/DECOUPLE/log/Mod1_%d_%A_%a.log
#SBATCH --cpus-per-task=3
#SBATCH --time=20:00:00
#SBATCH --mem=30000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=57-112
##SBATCH --array=1-112
##SBATCH --array=1-56

### %A is monsoon job number %a is interior array index

module load R/4.1.2 # load a specific R version

chmod +x shell_scripts/mod1_job.sh # for permissions
chmod +x scripts/02a_run_model1.R # for permissions

##site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteEND)
varname=$(sed -n "$SLURM_ARRAY_TASK_ID"p varnameEND)
sitename=$(sed -n "$SLURM_ARRAY_TASK_ID"p sitenameEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedEND)
scale=$(sed -n "$SLURM_ARRAY_TASK_ID"p scaleEND)

# Run the analysis
srun ./shell_scripts/mod1_job.sh $varname $sitename $seed $scale
