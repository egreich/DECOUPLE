#!/bin/bash
#SBATCH --job-name=DEPART_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/DECOUPLE/log/DEPART_%A_%a.log
#SBATCH --cpus-per-task=3
#SBATCH --time=8:00:00
#SBATCH --mem=40000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=1-9

### %A is monsoon job number %a is interior array index

module load R/4.1.2 # load a specific R version

chmod +x shell_scripts/DEPART_job.sh # for permissions
chmod +x scripts/02_run_DEPART.R # for permissions

##site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteEND)
sitename=$(sed -n "$SLURM_ARRAY_TASK_ID"p DEPARTsitenameEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p DEPARTseedEND)

# Run the analysis
srun ./shell_scripts/DEPART_job.sh $sitename $seed
