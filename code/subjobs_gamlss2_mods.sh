#!/bin/bash
#
#SBATCH --job-name=gamlss2
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --array=24,40
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/R-%A_%a.err

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts

# --- robust parse + debugging ---
CONFIGFN="$1"
echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

# -------------------------------------------------------------------
# Parse config line for this array task
line=$(awk -F'\t' -v id="$SLURM_ARRAY_TASK_ID" '$1==id {print; exit}' "$CONFIGFN" | tr -d '\r')

if [ -z "$line" ]; then
  echo "ERROR: no matching line for ID $SLURM_ARRAY_TASK_ID in $CONFIGFN" >&2
  exit 1
fi

# Split on tabs
IFS=$'\t' read -r idx PHENO FORM NAME <<< "$line"

# Remove surrounding quotes from FORM
FORM=$(printf '%s' "$FORM" | sed 's/^"\(.*\)"$/\1/')

# Trim whitespace
PHENO=$(printf '%s' "$PHENO" | xargs)
FORM=$(printf '%s' "$FORM" | xargs)
NAME=$(printf '%s' "$NAME" | xargs)

# Print for sanity check
printf 'IDX:   [%s]\n' "$idx"
printf 'PHENO: [%s]\n' "$PHENO"
printf 'FORM:  [%s]\n' "$FORM"
printf 'NAME:  [%s]\n' "$NAME"
# -------------------------------------------------------------------

SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.2.3.sif"
script="$BASE/code/test_gamlss2.R"

echo "SCRIPT: $script"

singularity run --cleanenv \
    -B "$BASE" \
    "$SINGULARITY_IMAGE" \
    Rscript "$script" "$PHENO" "$FORM" "$NAME"

echo "Job finished running!"