#!/bin/bash
#
#SBATCH --job-name=datasplit
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/datasplit_%A_%a.out
#SBACTH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/datasplit_%A_%a.err

CONFIGFN=$1

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $1}' $CONFIGFN )
PHENO_LIST=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
NAME_PREFIX=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )

#------------------

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/

module load R/4.2.3

Rscript $BASE/code/split_data_for_combat.R $DF $PHENO_LIST $SAVE_PATH $NAME_PREFIX 

# Done!
echo "Job finished running!"
