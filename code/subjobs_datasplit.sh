#!/bin/bash
#
#SBATCH --job-name=datasplit
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/combat/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/combat/R-%A_%a.err
#SBATCH --array=1-5

CONFIGFN=$1

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
PHENO_LIST=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
BATCH_VAR=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )

echo "DF: $DF"
echo "PHENO_LIST: $PHENO_LIST"
echo "SAVE_PATH: $SAVE_PATH"
echo "BATCH_VAR: $BATCH_VAR"

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts

module load R/4.4.0

Rscript $BASE/code/split_data_for_combat.R $DF $PHENO_LIST $SAVE_PATH $BATCH_VAR 

# Done!
echo "Job finished running!"
