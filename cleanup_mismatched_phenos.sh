#!/usr/bin/env bash
# For each train .rds path listed in pheno_mismatches.txt, remove that file
# and the corresponding test artifacts in the mirrored cv_sample_<opp>_test
# directory, so that the weighted (or age2plus) train+test pipelines can be
# re-run for just those phenotypes with --rerun TRUE.
#
# Dry-run by default. Pass --apply to actually delete.
#
# Usage:
#   ./cleanup_mismatched_phenos.sh [--apply] [--mode weighted|age2plus] \
#       [--mismatches pheno_mismatches.txt] [--pheno-lists pheno_lists] \
#       [--root .]
#
# Each input path is expected to look like:
#   <root>/cv_sample_<X>_train/<MODE>_<spec>/model_objs/<pheno>_<rest>_train_<MODE>.rds
#
# Removed (where they exist):
#   <root>/cv_sample_<opp>_test/<MODE>_<spec>/model_objs/<basetest>_full_mod.rds
#   <root>/cv_sample_<opp>_test/<MODE>_<spec>/model_objs/<basetest>_null_mod.rds
#   <root>/cv_sample_<opp>_test/<MODE>_<spec>/model_sums/<basetest>_summary.csv
#   <root>/cv_sample_<opp>_test/<MODE>_<spec>/model_sums/<basetest>_LRtest.csv
#   <root>/cv_sample_<opp>_test/<MODE>_<spec>/cent_csvs/<basetest>_centiles.csv
#   <root>/cv_sample_<opp>_test/<MODE>_<spec>/centile_plots/<basetest>.png
#   <root>/cv_sample_<opp>_test/<MODE>_<spec>/worm_plots/<basetest>.png
#   <root>/cv_sample_<opp>_test/<MODE>_<spec>/cent_csvs/<pheno>_sexdiffs.csv
#   <root>/cv_sample_<opp>_test/<MODE>_<spec>/cent_csvs/<pheno>_uniform_sexdiffs.csv
# where <basetest> = <pheno>_<rest>_test_<MODE>

set -u
shopt -s nullglob

mode="weighted"
apply=0
mismatches="pheno_mismatches.txt"
pheno_dir="pheno_lists"
root="."

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --apply) apply=1; shift ;;
        --mode) mode="${2:?}"; shift 2 ;;
        --mismatches) mismatches="${2:?}"; shift 2 ;;
        --pheno-lists) pheno_dir="${2:?}"; shift 2 ;;
        --root) root="${2:?}"; shift 2 ;;
        -h|--help) sed -n '1,30p' "$0"; exit 0 ;;
        *) echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done

case "$mode" in
    weighted|age2plus) ;;
    *) echo "unknown mode: $mode" >&2; exit 1 ;;
esac

file_suffix="_train_${mode}.rds"

if [[ ! -f "$mismatches" ]]; then
    echo "mismatches file not found: $mismatches" >&2
    exit 1
fi
if [[ ! -d "$pheno_dir" ]]; then
    echo "pheno list dir not found: $pheno_dir" >&2
    exit 1
fi

# Sort phenos by length desc so longer names match first (e.g. sGMV before GMV)
mapfile -t phenos < <(
    cat "${pheno_dir}"/*.txt \
        | sed '/^$/d' \
        | sort -u \
        | awk '{print length, $0}' | sort -k1,1nr -k2,2 | cut -d' ' -f2-
)
if (( ${#phenos[@]} == 0 )); then
    echo "no phenos found in ${pheno_dir}/*.txt" >&2
    exit 1
fi

extract_pheno() {
    local base="$1" p
    for p in "${phenos[@]}"; do
        if [[ "$base" == "${p}_"* ]]; then
            printf '%s\n' "$p"
            return 0
        fi
    done
    return 1
}

rm_path() {
    local p="$1"
    if [[ ! -e "$p" ]]; then
        printf '  miss   %s\n' "$p"
        return
    fi
    if (( apply )); then
        if rm -f -- "$p"; then
            printf '  rm     %s\n' "$p"
        else
            printf '  FAILED %s\n' "$p"
        fi
    else
        printf '  would-rm %s\n' "$p"
    fi
}

(( apply )) || echo "[dry-run] no files will be deleted. re-run with --apply to remove."
echo

total_lines=0
total_train=0
total_test_files=0
skipped=0

while IFS= read -r train_path || [[ -n "$train_path" ]]; do
    [[ -z "${train_path// /}" ]] && continue
    total_lines=$((total_lines + 1))

    base=$(basename "$train_path")
    if [[ "$base" != *"$file_suffix" ]]; then
        echo "skip (suffix mismatch for mode=$mode): $train_path" >&2
        skipped=$((skipped + 1))
        continue
    fi

    # Path components
    # .../cv_sample_<X>_train/<MODE>_<spec>/model_objs/<base>
    train_mo_dir=$(dirname "$train_path")           # .../model_objs
    spec_dir=$(dirname "$train_mo_dir")             # .../<MODE>_<spec>
    cv_dir=$(dirname "$spec_dir")                   # .../cv_sample_<X>_train
    cv_name=$(basename "$cv_dir")                   # cv_sample_<X>_train
    spec_name=$(basename "$spec_dir")               # <MODE>_<spec>

    # split letter (A/B/...) -> opposite
    if [[ "$cv_name" =~ ^cv_sample_(.+)_train$ ]]; then
        split="${BASH_REMATCH[1]}"
    else
        echo "skip (cannot parse cv_sample from): $train_path" >&2
        skipped=$((skipped + 1))
        continue
    fi
    case "$split" in
        A) opp="B" ;;
        B) opp="A" ;;
        *)  echo "skip (unexpected split '$split', not A/B): $train_path" >&2
            skipped=$((skipped + 1)); continue ;;
    esac

    # Extract pheno
    if ! pheno=$(extract_pheno "$base"); then
        echo "skip (no pheno prefix match): $train_path" >&2
        skipped=$((skipped + 1))
        continue
    fi

    # Build test-side names
    test_cv_dir="$(dirname "$cv_dir")/cv_sample_${opp}_test"
    test_spec_dir="${test_cv_dir}/${spec_name}"
    base_no_ext="${base%.rds}"                       # <pheno>_<rest>_train_<MODE>
    test_base="${base_no_ext/_train_/_test_}"        # <pheno>_<rest>_test_<MODE>

    echo "=== $train_path ==="
    echo "  pheno=$pheno  split=$split  opp=$opp  spec=$spec_name"

    # Train file
    echo "  [train]"
    rm_path "$train_path"
    total_train=$((total_train + 1))

    # Test files
    echo "  [test artifacts under ${test_spec_dir}]"
    declare -a victims=(
        "${test_spec_dir}/model_objs/${test_base}_full_mod.rds"
        "${test_spec_dir}/model_objs/${test_base}_null_mod.rds"
        "${test_spec_dir}/model_sums/${test_base}_summary.csv"
        "${test_spec_dir}/model_sums/${test_base}_LRtest.csv"
        "${test_spec_dir}/cent_csvs/${test_base}_centiles.csv"
        "${test_spec_dir}/centile_plots/${test_base}.png"
        "${test_spec_dir}/worm_plots/${test_base}.png"
        "${test_spec_dir}/cent_csvs/${pheno}_sexdiffs.csv"
        "${test_spec_dir}/cent_csvs/${pheno}_uniform_sexdiffs.csv"
    )
    for v in "${victims[@]}"; do
        rm_path "$v"
        total_test_files=$((total_test_files + 1))
    done
    echo
done < "$mismatches"

echo "----"
echo "lines processed: $total_lines (skipped: $skipped)"
echo "train .rds targeted: $total_train"
echo "test artifacts targeted: $total_test_files"
echo

if (( apply )); then
    echo "Done. To rerun the ${mode} train+test pipelines for just the deleted phenos,"
    echo "use --rerun TRUE so already-fit phenos are skipped:"
else
    echo "Re-run this script with --apply to actually delete."
    echo "After deletion, kick off the ${mode} train+test pipelines with --rerun TRUE:"
fi
echo
echo "  # adjust --total / --log_age to match your run"
if [[ "$mode" == "weighted" ]]; then
    cat <<'EOF'
  # 1. (re)build rerun config files for train and test
  bash code/config_cv_mods_weighted.sh      --total TRUE --log_age TRUE --rerun TRUE
  bash code/config_cv_mods_test_weighted.sh --total TRUE --log_age TRUE --rerun TRUE

  # 2. submit the SBATCH array jobs against the newly written rerun configs
  #    (config file path is printed by the config scripts; example pattern:)
  #    code/config_files/cv_sample_<A|B>_total<T>_logAge<L>_weighted_rerun<YYYYMMDD>_config.txt
  for split in A B; do
      sbatch code/subjobs_cv_mods_weighted.sh \
          code/config_files/cv_sample_${split}_totalTRUE_logAgeTRUE_weighted_rerun$(date +%Y%m%d)_config.txt
      sbatch code/subjobs_cv_mods_test_weighted.sh \
          code/config_files/cv_sample_${split}_totalTRUE_logAgeTRUE_weighted_test_rerun$(date +%Y%m%d)_config.txt
  done
EOF
else
    cat <<'EOF'
  # 1. (re)build rerun config files for train and test
  bash code/config_cv_mods_age2plus.sh      --total TRUE --log_age TRUE --rerun TRUE
  bash code/config_cv_mods_test_age2plus.sh --total TRUE --log_age TRUE --rerun TRUE

  # 2. submit the SBATCH array jobs against the newly written rerun configs
  for split in A B; do
      sbatch code/subjobs_cv_mods_age2plus.sh \
          code/config_files/cv_sample_${split}_totalTRUE_logAgeTRUE_age2plus_rerun$(date +%Y%m%d)_config.txt
      sbatch code/subjobs_cv_mods_test_age2plus.sh \
          code/config_files/cv_sample_${split}_totalTRUE_logAgeTRUE_age2plus_test_rerun$(date +%Y%m%d)_config.txt
  done
EOF
fi
