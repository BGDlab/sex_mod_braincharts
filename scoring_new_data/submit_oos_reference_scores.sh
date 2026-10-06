#!/bin/bash
# Score new data with oos_reference_scores.R as a SLURM job array (one phenotype per task),
# then stitch the per-phenotype outputs into one file with a dependent job.
#
# usage (run on the login node, NOT through sbatch):
#   bash scoring_new_data/submit_oos_reference_scores.sh \
#     --df my_datadict.csv --total TRUE --out_file my_ref_scores.csv \
#     [any other oos_reference_scores.R options except --pheno_list]
#
# every phenotype in all_phenos.txt that is a column of --df gets scored
#
# environment variables (optional):
#   RSCRIPT       command used to run R, e.g.
#                 RSCRIPT="singularity run --cleanenv -B /mnt/isilon img.sif Rscript"
#   MAX_PARALLEL  max array tasks running at once (default 50)
#   SBATCH_ARGS   extra sbatch flags for the array tasks, e.g. "--mem=16G --time=2:00:00"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"
MAX_PARALLEL="${MAX_PARALLEL:-50}"
SBATCH_ARGS="${SBATCH_ARGS:-}"

# pull --out_file out of the args (and note --df); pass everything else through to R
OUT_FILE=""
DF=""
PASS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out_file|--out-file)       OUT_FILE="$2"; shift 2 ;;
    --out_file=*|--out-file=*)   OUT_FILE="${1#*=}"; shift ;;
    --pheno_list*|--pheno-list*)
      echo "--pheno_list isn't used here: every phenotype in all_phenos.txt found in --df is scored" >&2
      exit 1 ;;
    --df)   DF="$2"; PASS+=("$1" "$2"); shift 2 ;;
    --df=*) DF="${1#*=}"; PASS+=("$1"); shift ;;
    *) PASS+=("$1"); shift ;;
  esac
done
[[ -z "$OUT_FILE" ]] && { echo "missing required option: --out_file" >&2; exit 1; }
[[ -z "$DF" ]] && { echo "missing required option: --df" >&2; exit 1; }

OUT_FILE="$(cd "$(dirname "$OUT_FILE")" && pwd)/$(basename "$OUT_FILE")"
WORK_DIR="${OUT_FILE%.*}_array"
mkdir -p "$WORK_DIR"/{scores,pheno_lists,logs}

# full list of phenotypes with models: local copy next to this script, else from GitHub
ALL_PHENOS="$SCRIPT_DIR/all_phenos.txt"
if [[ ! -f "$ALL_PHENOS" ]]; then
  ALL_PHENOS="$WORK_DIR/all_phenos.txt"
  curl -fsSL "https://raw.githubusercontent.com/BGDlab/sex_mod_braincharts/refs/heads/sharing/scoring_new_data/all_phenos.txt" \
    -o "$ALL_PHENOS"
fi

# keep the phenotypes that are columns of --df, so every task (and the combine step)
# sees the same list
# shellcheck disable=SC2086
$RSCRIPT -e '
  a <- commandArgs(trailingOnly = TRUE)
  phenos <- trimws(readLines(a[1], warn = FALSE))
  phenos <- phenos[nzchar(phenos)]
  cols <- names(data.table::fread(a[2], nrows = 0))
  keep <- intersect(phenos, cols)
  writeLines(keep, a[3])
  cat(length(keep), "/", length(phenos), " phenotypes found in ", a[2], "\n", sep = "")
' "$ALL_PHENOS" "$DF" "$WORK_DIR/phenos.txt"
N=$(grep -c . "$WORK_DIR/phenos.txt" || true)
[[ "$N" -eq 0 ]] && { echo "none of the phenotypes in all_phenos.txt are columns of $DF" >&2; exit 1; }

# install/update packages once here, so the array tasks don't all try to at once
# shellcheck disable=SC2086
$RSCRIPT "$SCRIPT_DIR/oos_reference_scores.R" --install_only

export OOS_SCRIPT_DIR="$SCRIPT_DIR" OOS_WORK_DIR="$WORK_DIR" RSCRIPT

# shellcheck disable=SC2086
ARRAY_ID=$(sbatch --parsable \
  --array=1-"$N"%"$MAX_PARALLEL" \
  --output="$WORK_DIR/logs/score_%A_%a.out" \
  $SBATCH_ARGS \
  "$SCRIPT_DIR/subjob_oos_reference_scores.sh" "${PASS[@]}")

# afterany: still combine if some tasks fail, so the missing phenos get reported
COMBINE_ID=$(sbatch --parsable \
  --job-name=oos_combine \
  --dependency=afterany:"$ARRAY_ID" \
  --time=1:00:00 --mem=16G \
  --output="$WORK_DIR/logs/combine_%j.out" \
  --wrap="$RSCRIPT '$SCRIPT_DIR/combine_oos_reference_scores.R' '$WORK_DIR' '$OUT_FILE'")

echo "submitted array job $ARRAY_ID ($N phenotypes) and combine job $COMBINE_ID"
echo "per-phenotype outputs & logs: $WORK_DIR"
echo "final output: $OUT_FILE"
