#!/bin/bash
#
#SBATCH --job-name=check_converge
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-20

CONFIGFN=$1
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
echo "DF: $DF"

# Singularity image path 
SINGULARITY_IMAGE="/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/containers/r_gamlss_0.2.4.sif"

#------------------
BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/

# Run R script inside Singularity container
singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $BASE/code/check_convergence.R $DF

# Done!
echo "Job finished running!"
