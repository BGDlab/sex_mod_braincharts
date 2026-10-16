#!/bin/bash
# Score new data with oos_reference_scores.R as a SLURM job array (one phenotype per task),
# then stitch the per-phenotype outputs into one file with a dependent job.
#
# usage:
#   bash scoring_new_data/submit_oos_reference_scores.sh \
#     --df my_data.csv --total TRUE --out_file my_ref_scores.csv \
#     [any other oos_reference_scores.R options except --pheno_list]
#
# every phenotype in all_phenos.txt that is a column of --df gets scored
#
# environment variables (optional):
#   MAX_PARALLEL  max array tasks running at once (default 50)
#   SBATCH_ARGS   extra sbatch flags for the array tasks, e.g. "--mem=16G --time=2:00:00"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAX_PARALLEL="${MAX_PARALLEL:-50}"
SBATCH_ARGS="${SBATCH_ARGS:-}"

# pull --out_file out of the args (and note --df and --ref_data); pass everything else through to R
OUT_FILE=""
DF=""
REF_DATA=""
PASS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out_file|--out-file)       OUT_FILE="$2"; shift 2 ;;
    --out_file=*|--out-file=*)   OUT_FILE="${1#*=}"; shift ;;
    --pheno_list*|--pheno-list*)
      echo "--pheno_list isn't used here: every phenotype in all_phenos.txt found in --df is scored" >&2
      exit 1 ;;
    --df)   DF="$2"; shift 2 ;;
    --df=*) DF="${1#*=}"; shift ;;
    --ref_data|--ref-data)     REF_DATA="$2"; PASS+=("$1" "$2"); shift 2 ;;
    --ref_data=*|--ref-data=*) REF_DATA="${1#*=}"; PASS+=("$1"); shift ;;
    *) PASS+=("$1"); shift ;;
  esac
done
[[ -z "$OUT_FILE" ]] && { echo "missing required option: --out_file" >&2; exit 1; }
[[ -z "$DF" ]] && { echo "missing required option: --df" >&2; exit 1; }

[[ -d "$(dirname "$OUT_FILE")" ]] || { echo "--out_file folder not found: $(dirname "$OUT_FILE")" >&2; exit 1; }
OUT_FILE="$(cd "$(dirname "$OUT_FILE")" && pwd)/$(basename "$OUT_FILE")"
WORK_DIR="${OUT_FILE%.*}_array"
mkdir -p "$WORK_DIR"/{scores,logs}

# every task and the combine step read a frozen copy of --df (written below) with the rows
# numbered in advance (.row_id), so editing the original while jobs are queued can't misalign rows
[[ "$DF" =~ ^https?:// || -f "$DF" ]] || { echo "--df file not found: $DF" >&2; exit 1; }
DF_COPY="$WORK_DIR/df.csv"
PASS+=(--df "$DF_COPY")

# full list of phenotypes with models: local copy next to this script, else from GitHub
ALL_PHENOS="$SCRIPT_DIR/all_phenos.txt"
if [[ ! -f "$ALL_PHENOS" ]]; then
  ALL_PHENOS="$WORK_DIR/all_phenos.txt"
  curl -fsSL "https://raw.githubusercontent.com/BGDlab/sex_mod_braincharts/refs/heads/sharing/scoring_new_data/all_phenos.txt" \
    -o "$ALL_PHENOS"
fi

# check --df has the columns every model needs (and any named in --ref_data), write the
# numbered copy, then list out scorable phenos (i.e. all phenos in --df)
Rscript -e '
  a <- commandArgs(trailingOnly = TRUE)
  source(a[4])
  df <- data.table::fread(a[2], na.strings = c("NA", "", "\"\""))
  check_df(df, if (nzchar(a[5])) a[5])
  if (".row_id" %in% names(df))
    stop("--df already has a .row_id column; rename or remove it", call. = FALSE)
  data.table::set(df, j = ".row_id", value = seq_len(nrow(df)))
  data.table::setcolorder(df, ".row_id")
  data.table::fwrite(df, a[6])
  phenos <- trimws(readLines(a[1], warn = FALSE))
  phenos <- phenos[nzchar(phenos)]
  keep <- intersect(phenos, names(df))
  writeLines(keep, a[3])
  cat(length(keep), "/", length(phenos), " phenotypes found in ", a[2], "\n", sep = "")
' "$ALL_PHENOS" "$DF" "$WORK_DIR/phenos.txt" "$SCRIPT_DIR/oos_reference_scores_helper_funs.R" "$REF_DATA" "$DF_COPY"
N=$(grep -c . "$WORK_DIR/phenos.txt" || true)
[[ "$N" -eq 0 ]] && { echo "none of the phenotypes in all_phenos.txt are columns of $DF" >&2; exit 1; }

# install/update packages once
Rscript "$SCRIPT_DIR/oos_reference_scores.R" --install_only

export OOS_SCRIPT_DIR="$SCRIPT_DIR" OOS_WORK_DIR="$WORK_DIR"

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
  --wrap="Rscript '$SCRIPT_DIR/combine_oos_reference_scores.R' '$WORK_DIR' '$OUT_FILE'")

echo "submitted array job $ARRAY_ID ($N phenotypes) and pending combine job $COMBINE_ID"
echo "per-phenotype outputs & logs: $WORK_DIR"
echo "final output: $OUT_FILE"
