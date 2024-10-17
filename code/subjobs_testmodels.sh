#!/bin/bash
#
#SBATCH --job-name=testmodels
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-224
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/testmodels_%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/testmodels_%A_%a.err

CONFIGFN=$1

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
PHENO=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
KNOTS=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )

echo "DF: $DF"
echo "PHENO: $PHENO"
echo "KNOTS: $KNOTS"
echo "SAVE_PATH: $SAVE_PATH"

#------------------

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/

module load R/4.4.0

Rscript $BASE/code/fit_test_models.R $DF $PHENO $KNOTS $SAVE_PATH 

# Done!
echo "Job finished running!"
