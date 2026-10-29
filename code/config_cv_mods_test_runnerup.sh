#!/bin/bash
#script to write config file that will test runner-up (2nd-lowest BIC) training models on other sample
#parallel check on config_cv_mods_test.sh - same config format, so run w/ subjobs_cv_mods_test.sh
#requires manifest from select_runnerup_mods.R
#run from outside code dir

#PATHS
data_path=./data
config_path=./code/config_files
#hardcoding sim df for consistency across phenos and splits
sim_df=$(realpath $data_path/v3_full_sim.rds)

rerun="FALSE"

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --total) total="$2"; shift 2 ;;
    --log_age) log_age="$2"; shift 2 ;;
    --rerun) rerun="$2"; shift 2 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

# Check that all required arguments are provided
if [[ -z "$total" || -z "$log_age" ]]; then
  echo "Missing arguments. Usage:"
  echo "--total TRUE/FALSE --log_age TRUE/FALSE [--rerun TRUE/FALSE]"
  exit 1
fi

# Print arguments
echo "include total value = $total"
echo "log scale age = $log_age"

manifest=$config_path/cv_runnerup_manifest_total${total}_logAge${log_age}.csv
if ! [ -f $manifest ]; then
  echo "Manifest not found: $manifest"
  echo "Run first: Rscript code/select_runnerup_mods.R $total $log_age"
  exit 1
fi

#LOOP THROUGH 1/2 CSVS (split = sample being tested in)
for split in A B
do

  echo "prepping: $split"

    #make config file dir or remove old file if necessary
    if [[ "$rerun" == "TRUE" ]]; then
      date_tag=$(date +%Y%m%d)
      config_file=$config_path/cv_sample_${split}_total${total}_logAge${log_age}_test_runnerup_rerun${date_tag}_config.txt
    else
      config_file=$config_path/cv_sample_${split}_total${total}_logAge${log_age}_test_runnerup_config.txt
    fi
    if ! [ -d $config_path ]
    then
      mkdir $config_path
    elif [ -f $config_file ]
    then
      rm -rf $config_file
    fi

    touch $config_file

  #make output dir - kept separate from cv_sample_?_test so downstream globs don't pick up runner-ups
  save_dir=./cv_sample_${split}_test_runnerup
  if ! [ -d $save_dir ]
  then
    mkdir $save_dir
  fi

  #LOOP THROUGH MANIFEST ROWS (skip header)
  while IFS=, read -r train_split test_split pheno_cat pheno status og_mod _rest
  do
    [[ "$test_split" == "$split" ]] || continue

    if [[ "$status" != "ok" ]]; then
      echo "Warning: skipping $pheno ($pheno_cat, trained on $train_split) - manifest status: $status" >&2
      continue
    fi

    #CREATE OUTPUT DIRS
    save_path=$save_dir/${pheno_cat}_total${total}_logAge${log_age}_pbmods
    if ! [ -d $save_path ]
    then
      mkdir $save_path
      mkdir $save_path/model_objs
      mkdir $save_path/centile_plots
      mkdir $save_path/worm_plots
      mkdir $save_path/cent_csvs
      mkdir $save_path/model_sums
    fi
    save_path=$(realpath $save_path)

    #skip already tested models if rerunning
    if [[ "$rerun" == "TRUE" && -f ${save_path}/cent_csvs/${pheno}_uniform_sexdiffs.csv ]]; then
      echo "Skipping $pheno (uniform sex diffs csv found)"
      continue
    fi

    #write csv to test in - handle optional _logPheno*_ in filename
    mapfile -t file_matches < <(find $(realpath $data_path/cv_sample_${split}_dfs) -type f -name "${pheno}_total${total}*logAge${log_age}.csv" 2>/dev/null)
    if [ ${#file_matches[@]} -gt 1 ]; then
      echo "Error: Multiple CSV files found for $pheno:"
      printf '%s\n' "${file_matches[@]}"
      exit 1
    elif [ ${#file_matches[@]} -eq 0 ]; then
      echo "Warning: No CSV found for $pheno, skipping"
      continue
    else
      file="${file_matches[0]}"
    fi

    # Write the CSV file path and the model path to the output file (tab-delimited)
    echo -e "$file\t$og_mod\t$save_path\t$total\t$sim_df" >> "$config_file"
  done < <(tail -n +2 "$manifest")

  #add numbering
  nl "$config_file" > temp.txt && mv temp.txt "$config_file"
  echo "$(wc -l < $config_file) models to test in $config_file"
done
