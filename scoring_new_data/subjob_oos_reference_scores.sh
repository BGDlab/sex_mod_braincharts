#!/bin/bash
# One task of the array submitted by submit_oos_reference_scores.sh: scores a single phenotype.
# Args are passed straight through to oos_reference_scores.R.
#
#SBATCH --job-name=oos_scores
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

PHENO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$OOS_WORK_DIR/phenos.txt")
echo "TASK $SLURM_ARRAY_TASK_ID: $PHENO"

PHENO_FILE="$OOS_WORK_DIR/pheno_lists/$PHENO.txt"
echo "$PHENO" > "$PHENO_FILE"

# shellcheck disable=SC2086
$RSCRIPT "$OOS_SCRIPT_DIR/oos_reference_scores.R" "$@" \
  --skip_install \
  --pheno_list "$PHENO_FILE" \
  --out_file "$OOS_WORK_DIR/scores/$PHENO.csv"

echo "Job finished running!"
