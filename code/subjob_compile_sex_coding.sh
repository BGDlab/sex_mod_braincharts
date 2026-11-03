#!/bin/bash
#
#SBATCH --job-name=sex_coding_pdf
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/code/jobfiles/sex_coding/compile-%j.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/code/jobfiles/sex_coding/compile-%j.err

#compile refit_sex_coding.R outputs for every BestMod in the list into one pdf
BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision

LIST=$1
OUT_PDF=$2

echo "Model list: $LIST"
echo "PDF: $OUT_PDF"

SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.2.14.sif"

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $BASE/code/compile_sex_coding_pdf.R $LIST $OUT_PDF

echo "Job finished running!"
