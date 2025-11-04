#!/bin/bash
#
#SBATCH --job-name=replots
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

# Singularity image path 
SINGULARITY_IMAGE="/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/containers/r_gamlss_0.2.4.sif"

#------------------
BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/

# Load Singularity module (adjust module name as needed for your cluster)
#module load singularity

# Run R script inside Singularity container
# Bind mount the necessary directories so the container can access your data
singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $BASE/code/check_convergence.R

# Done!
echo "Job finished running!"
