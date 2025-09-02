#!/bin/bash
#
#SBATCH --job-name=gamlss2
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-8
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/R-%A_%a.err

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts
CONFIGFN=$1

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
PHENO=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' "$CONFIGFN")
FORM=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' "$CONFIGFN")
RAND=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' "$CONFIGFN")

echo "PHENO: $PHENO"
echo "FORM: $FORM"
echo "RAND: $RAND"

#------------------

SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.2.2.sif"

script=$BASE/code/test_gamlss2.R

echo "SCRIPT: $script"

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $script "$PHENO" "$FORM" "$RAND"

# Done!
echo "Job finished running!"
