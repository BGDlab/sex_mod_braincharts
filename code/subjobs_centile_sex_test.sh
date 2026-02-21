#!/bin/bash
#
#SBATCH --job-name=cent_sex_test
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15G
#SBATCH --array=1-1446%50
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_test/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_test/R-%A_%a.err

CONFIGFN=$1
#CONFIGFN=$(realpath $CONFIGFN)

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
MODEL=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
DX=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )

echo "DF: $DF"
echo "TRAINING MODEL: $MODEL"
echo "SAVE_PATH: $SAVE_PATH"
echo "DX: $DX"

#------------------

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/
SINGULARITY_IMAGE="/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/containers/r_gamlss_0.2.11.sif"

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $BASE/code/centile_sex_test.R $DF $MODEL $SAVE_PATH $DX

# Done!
echo "Job finished running!"
