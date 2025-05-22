#!/bin/bash
#
#SBATCH --job-name=replots
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_train/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/cv_train/R-%A_%a.err

# Singularity image path 
SINGULARITY_IMAGE="/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/containers/r_gamlss_0.1.0.sif"

#------------------
BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/

# Load Singularity module (adjust module name as needed for your cluster)
#module load singularity

# Run R script inside Singularity container
# Bind mount the necessary directories so the container can access your data
singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $BASE/code/replot_centiles.R

# Done!
echo "Job finished running!"
