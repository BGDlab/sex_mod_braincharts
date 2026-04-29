#!/bin/bash

PREFIX_FILE="sczglu_pheno_rerun.txt"

# Read prefixes into an array (strips any trailing whitespace/carriage returns)
mapfile -t prefixes < <(tr -d '\r' < "$PREFIX_FILE")

found=0

for dir in cv_sample_?_test/*/*/; do
    [ -d "$dir" ] || continue

    subdir=$(basename "$(dirname "$dir")")

    # Only process subdirs starting with the three target strings
    if [[ "$subdir" != cortical_vols* && "$subdir" != global_vols* && "$subdir" != subcortical* ]]; then
        continue
    fi

    for file in "$dir"*; do
        [ -f "$file" ] || continue
        fname=$(basename "$file")

        for prefix in "${prefixes[@]}"; do
            if [[ "$fname" == "$prefix"* ]]; then
                echo "remove $file"
		rm "$file"
                ((found++))
                break
            fi
        done
    done
done

echo ""
echo "Total files that would be removed: $found"
