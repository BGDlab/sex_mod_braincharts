#!/bin/bash
#
#SBATCH --job-name=cv_mods
#SBATCH --time=336:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-4%10
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_train/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_train/R-%A_%a.err

CONFIGFN=$1

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
PHENO=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
FS=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
TOTAL=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $CONFIGFN )
LOG_PHENO=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $CONFIGFN )
LOG_AGE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $CONFIGFN )
SMOOTH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $9}' $CONFIGFN )

echo "DF: $DF"
echo "PHENO: $PHENO"
echo "FREESURFER COV: $FS"
echo "TOTAL COV: $TOTAL"
echo "SAVE_PATH: $SAVE_PATH"
echo "LOG-SCALE PHENO: $LOG_PHENO"
echo "LOG-SCALE AGE: $LOG_AGE"
echo "SMOOTH: $SMOOTH"

#------------------

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/

if [ $TOTAL == "NULL" ] 
then
  script=$BASE/code/fit_cv_mods.R
else
  script=$BASE/code/fit_cv_total_mods.R
fi

echo "SCRIPT: $script"

module load R/4.4.0

Rscript $script $DF $PHENO $LAMBDA $FS $TOTAL $SAVE_PATH $WEIGHT $LOG_PHENO $LOG_AGE $SMOOTH

# Done!
echo "Job finished running!"
