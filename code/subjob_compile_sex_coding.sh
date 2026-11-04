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

#run outside the gamlss container (no png package there), same as subjob_grab_pngs.sh
module load R/4.4.0
Rscript $BASE/code/compile_sex_coding_pdf.R $LIST $OUT_PDF

echo "Job finished running!"
