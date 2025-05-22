#!/bin/bash
#
#SBATCH --job-name=cv_dfs
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-224%20 #total count = 224
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_dfs/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_dfs/R-%A_%a.err

CONFIGFN=$1
#CONFIGFN=$(realpath $CONFIGFN)

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
PHENO=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
FS=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
TOTAL=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )
LOG_PHENO=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $CONFIGFN )
LOG_AGE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $CONFIGFN )

echo "DF: $DF"
echo "PHENOTYPE: $PHENO"
echo "FREESURFER: $FS"
echo "TOTAL: $TOTAL"
echo "SCALE PHENO: $LOG_PHENO"
echo "SCALE AGE: $LOG_AGE"
echo "SAVE_PATH: $SAVE_PATH"

#------------------

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/

module load R/4.4.0

Rscript $BASE/code/prep_cv_dfs.R $DF $PHENO $FS $TOTAL $LOG_PHENO $LOG_AGE $SAVE_PATH

# Done!
echo "Job finished running!"
