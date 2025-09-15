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
CONFIGFN=$1

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

# --- read the first matching line (preserve everything) ---
line=$(awk -F'\t' -v id="$SLURM_ARRAY_TASK_ID" '$1==id {print; exit}' "$CONFIGFN")
# remove only Windows CR if present, but do NOT touch backslashes or quotes
line=$(printf '%s' "$line" | tr -d '\r')

if [ -z "$line" ]; then
  echo "ERROR: no matching line for ID $SLURM_ARRAY_TASK_ID in $CONFIGFN" >&2
  exit 1
fi

# split on TAB; -r prevents read from treating backslashes as escapes
IFS=$'\t' read -r idx PHENO FORM NAME <<< "$line"

# trim leading/trailing whitespace **without** altering backslashes/quotes
PHENO=$(printf '%s' "$PHENO" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
FORM=$(printf '%s' "$FORM" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
NAME=$(printf '%s' "$NAME" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')

# Diagnostics
printf 'IDX:   [%s]\n' "$idx"
printf 'PHENO: [%s]\n' "$PHENO"
printf 'FORM:  [%s]\n' "$FORM"
printf 'NAME:  [%s]\n' "$NAME"

# show raw bytes for FORM if you need to debug further:
printf '\nFORM BYTES (od -c):\n'; printf '%s' "$FORM" | od -c

# -------------------------------------------------------------------
SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.2.3.sif"
script="$BASE/code/test_gamlss2.R"

echo "SCRIPT: $script"

singularity run --cleanenv \
    -B "$BASE" \
    "$SINGULARITY_IMAGE" \
    Rscript "$script" "$PHENO" "$FORM" "$NAME"

echo "Job finished running!"
