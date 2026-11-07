#!/bin/bash
#
#SBATCH --job-name=sex_coding
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/code/jobfiles/sex_coding/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/code/jobfiles/sex_coding/R-%A_%a.err

#refit one BestMod (line SLURM_ARRAY_TASK_ID of the list) with sexFemale, sexEffect & pb + pvc(by = sexEffect) coding
#+ plot all 4 codings
#list: one BestMod path per line, optionally followed by a tab and its training csv
#submit with submit_sex_coding.sh

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision

LIST=$1

echo "Model list: $LIST"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE LIST (blank & # lines skipped, same as compile_sex_coding_pdf.R)
LINE=$(grep -v -e '^[[:space:]]*$' -e '^#' $LIST | sed -n "${SLURM_ARRAY_TASK_ID}p")
MODEL=$(echo "$LINE" | cut -f1)
CSV=$(echo "$LINE" | cut -s -f2)
CSV=${CSV:-NA}

echo "MODEL: $MODEL"
echo "CSV: $CSV"

#------------------

SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.2.14.sif"

script=$BASE/code/refit_sex_coding.R

echo "SCRIPT: $script"

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $script $MODEL $CSV

chmod 777 -R $(dirname $(dirname $MODEL))/sex_coding

# Done!
echo "Job finished running!"
