#!/bin/bash
#
#SBATCH --job-name=cv_mods
#SBATCH --time=80:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/code/jobfiles/cv_test/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/code/jobfiles/cv_test/R-%A_%a.err

CONFIGFN=$1
#CONFIGFN=$(realpath $CONFIGFN)

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
MODEL=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
TBV=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )
SIMDF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $CONFIGFN )

echo "DF: $DF"
echo "TRAINING MODEL: $MODEL"
echo "SAVE_PATH: $SAVE_PATH"
echo "TBV-corrected: $TBV"
echo "SIM DF: $SIMDF"

#------------------

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/
SINGULARITY_IMAGE="/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/containers/r_gamlss_0.2.9.sif"
export SINGULARITYENV_R_LIBS=$BASE/Rlib_extra:/usr/local/lib/R/site-library:/usr/local/lib/R/library #extra R pkgs not in image (ggnewscale); keep image libs ahead of ~/R

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $BASE/code/fit_cv_mods_test.R $DF $MODEL $SAVE_PATH $TBV $SIMDF

# Done!
echo "Job finished running!"
