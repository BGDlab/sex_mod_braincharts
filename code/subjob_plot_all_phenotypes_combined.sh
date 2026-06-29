#!/bin/bash
#
#SBATCH --job-name=combined_pdf
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/combined_pdf/R-%A.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/combined_pdf/R-%A.err

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts

SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.2.15.sif"

script=$BASE/code/plot_all_phenotypes_combined.R

echo "SCRIPT: $script"

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $script

# Done!
echo "Job finished running!"
