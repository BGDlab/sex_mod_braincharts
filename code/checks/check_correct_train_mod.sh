#!/usr/bin/env bash
# For every pheno listed in pheno_lists/*.txt, verify that each file
#   <cv>/<MODE>_<X>/model_objs/<pheno>_<suffix>_train_<MODE>.rds
# has a matching BestMod counterpart
#   <cv>/<X>/model_objs/<pheno>_<suffix>_BestMod.rds
# i.e. the part between "<pheno>_" and the trailing suffix is identical.
#
# MODE is either "weighted" (default) or "age2plus".
#
# Usage:
#   ./check_pheno_dupes.sh [--mode weighted|age2plus] \
#       [pheno_lists_dir] [root_dir] [out_txt]
#
# Defaults: --mode weighted  pheno_lists  .  pheno_mismatches.txt
#
# The .txt file lists each mismatched .rds path, one per line (no headers).
# The terminal still shows the detailed per-cv_sample report.

set -u
shopt -s nullglob

mode="weighted"
# parse optional --mode flag
if [[ "${1:-}" == "--mode" ]]; then
    mode="${2:?--mode requires a value (weighted|age2plus)}"
    shift 2
fi

case "$mode" in
    weighted|age2plus) ;;
    *) echo "unknown mode: $mode (expected weighted or age2plus)" >&2; exit 1 ;;
esac

pheno_dir="${1:-pheno_lists}"
root="${2:-.}"
out_txt="${3:-pheno_mismatches.txt}"

dir_prefix="${mode}_"           # e.g. weighted_  or  age2plus_
file_suffix="_train_${mode}.rds"  # e.g. _train_weighted.rds or _train_age2plus.rds

: > "$out_txt"

if [[ ! -d "$pheno_dir" ]]; then
    echo "pheno list dir not found: $pheno_dir" >&2
    exit 1
fi

mapfile -t phenos < <(cat "${pheno_dir}"/*.txt | sed '/^$/d' | sort -u)
if (( ${#phenos[@]} == 0 )); then
    echo "no phenos found in ${pheno_dir}/*.txt" >&2
    exit 1
fi

cv_dirs=( "${root}"/cv_sample_?_train )
if (( ${#cv_dirs[@]} == 0 )); then
    echo "no cv_sample_?_train dirs found under ${root}" >&2
    exit 1
fi

total_checked=0
total_mismatched=0

for cv in "${cv_dirs[@]}"; do
    [[ -d "$cv" ]] || continue
    cv_name=$(basename "$cv")
    echo "========================================"
    echo "=== $cv_name  (mode: $mode) ==="
    echo "========================================"
    cv_checked=0
    cv_mismatched=0

    for sdir in "$cv"/${dir_prefix}*/; do
        [[ -d "$sdir" ]] || continue
        sname=$(basename "$sdir")
        other_name="${sname#${dir_prefix}}"
        other_mo="${cv}/${other_name}/model_objs"

        for pheno in "${phenos[@]}"; do
            for sf in "${sdir}model_objs/${pheno}"_*"${file_suffix}"; do
                [[ -e "$sf" ]] || continue
                cv_checked=$((cv_checked + 1))
                sf_base=$(basename "$sf")
                expected_base="${sf_base%${file_suffix}}_BestMod.rds"
                expected_path="${other_mo}/${expected_base}"

                if [[ ! -f "$expected_path" ]]; then
                    cv_mismatched=$((cv_mismatched + 1))
                    echo "$sf" >> "$out_txt"
                    echo "MISMATCH pheno=$pheno"
                    echo "  ${mode}  : $sf"
                    echo "  expected : $expected_path"
                    existing=( "${other_mo}/${pheno}"_*_BestMod.rds )
                    if (( ${#existing[@]} > 0 )); then
                        echo "  found instead:"
                        for e in "${existing[@]}"; do echo "    $e"; done
                    else
                        echo "  (no BestMod files for this pheno in ${other_mo})"
                    fi
                fi
            done
        done
    done

    echo
    echo "-- $cv_name summary: checked $cv_checked ${mode} files, $cv_mismatched mismatches --"
    if (( cv_mismatched == 0 && cv_checked > 0 )); then
        echo "   all matched."
    fi
    echo

    total_checked=$((total_checked + cv_checked))
    total_mismatched=$((total_mismatched + cv_mismatched))
done

echo "========================================"
echo "overall (mode=$mode): checked $total_checked ${mode} files, $total_mismatched mismatches"
echo "mismatched paths written to: $out_txt"
(( total_mismatched == 0 )) && echo "All ${mode} files have matching BestMod counterparts."
