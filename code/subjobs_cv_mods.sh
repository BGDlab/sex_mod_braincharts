'#!/bin/bash
#
#SBATCH --job-name=cv_mods
#SBATCH --time=336:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_train/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_train/R-%A_%a.err

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts

CONFIGFN=$BASE/code/config_files/cv_mods_global_vols_logPhenoTRUE_totalTRUE_logAgeTRUE_smpb_config.txt

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

SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.1.2.sif"

script=$BASE/code/fit_cv_mods.R

echo "SCRIPT: $script"

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $script $DF $PHENO $FS $TOTAL $SAVE_PATH $LOG_PHENO $LOG_AGE $SMOOTH

chmod 777 -R $SAVE_PATH

# Done!
echo "Job finished running!"
