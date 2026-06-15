#!/bin/bash
#
#SBATCH --job-name=w_test
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1-241%100
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/w_test/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/w_test/R-%A_%a.err

CONFIGFN=$1

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
MODEL=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
TBV=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )
DF_SAVE_DIR=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $CONFIGFN )


echo "DF: $DF"
echo "TRAINING MODEL: $MODEL"
echo "SAVE_PATH: $SAVE_PATH"
echo "TBV-corrected: $TBV"
echo "DF SAVE DIR: $DF_SAVE_DIR"

#------------------

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/
SINGULARITY_IMAGE="/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/containers/r_gamlss_0.2.14.sif"

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $BASE/code/fit_cv_mods_test_weighted.R $DF $MODEL $SAVE_PATH $TBV $DF_SAVE_DIR

# Done!
echo "Job finished running!"
