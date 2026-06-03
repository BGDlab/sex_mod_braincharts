#!/bin/bash
#script to write a config file driving subjobs_replot.sh
#  one row per model = one slurm array task
#  writes one config per split (A, B)
#run from outside code dir

#PATHS
data_path=./data
config_path=./code/config_files
pheno_lists=./pheno_lists

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --total)      total="$2";      shift 2 ;;
    --log_age)    log_age="$2";    shift 2 ;;
    --train_test) train_test="$2"; shift 2 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

# Check that all required arguments are provided
if [[ -z "$total" || -z "$log_age" || -z "$train_test" ]]; then
  echo "Missing arguments. Usage:"
  echo "--total TRUE/FALSE --log_age TRUE/FALSE --train_test TRAIN/TEST"
  exit 1
fi

# Normalize: TRAIN/TEST (uppercase) is the label written to the config and used
# in plot filenames; train/test (lowercase) is what the project dirs are named.
train_test=$(echo "$train_test" | tr '[:lower:]' '[:upper:]')
tt_lower=$(echo "$train_test"   | tr '[:upper:]' '[:lower:]')

# Print arguments
echo "include total value = $total"
echo "log scale age = $log_age"
echo "train or test model = $train_test"

#LOOP THROUGH SPLITS
for split in A B
do

  echo "prepping: $split"

  #make config file dir or remove old file if necessary
  config_file=$config_path/cv_sample_${split}_total${total}_logAge${log_age}_${train_test}_replot_config.txt

  if ! [ -d $config_path ]
  then
    mkdir $config_path
  elif [ -f $config_file ]
  then
    rm -rf $config_file
  fi

  touch $config_file

  #LOOP THROUGH PHENO CATEGORIES
  for pheno_list in $(find $(realpath $pheno_lists) -type f -name "*.txt")
  do
    echo "pheno list: $pheno_list"

    #get filename -> freesurfer variable name
    pheno_cat=$(basename -- "$pheno_list")
    pheno_cat="${pheno_cat%.*}"

    #MODEL DIR (search-only; skip if it doesn't exist)
    save_dir=./cv_sample_${split}_${tt_lower}/${pheno_cat}_total${total}_logAge${log_age}_pbmods
    if ! [ -d $save_dir ]
    then
      echo "skipping missing $save_dir"
      continue
    fi
    save_dir=$(realpath $save_dir)

    #PLOT OUTPUT DIR (create if needed)
    save_path=$save_dir/replot
    if ! [ -d $save_path ]
    then
      mkdir -p $save_path
    fi
    save_path=$(realpath $save_path)

    #LOOP THROUGH PHENOS
    while read -r pheno_line
    do

      #find csv - handle optional _logPheno*_ in filename
      mapfile -t file_matches < <(find $(realpath $data_path/cv_sample_${split}_dfs) -type f -name "${pheno_line}_total${total}*logAge${log_age}.csv" 2>/dev/null)
      if [ ${#file_matches[@]} -gt 1 ]; then
        echo "Error: Multiple CSV files found for $pheno_line:"
        printf '%s\n' "${file_matches[@]}"
        exit 1
      elif [ ${#file_matches[@]} -eq 0 ]; then
        echo "Warning: No CSV found for $pheno_line, skipping"
        continue
      else
        file="${file_matches[0]}"
      fi

      #if train model, get BestModel
      if [[ "$train_test" == "TRAIN" ]]; then
        mapfile -t matches < <(find "$save_dir" -path "*model_objs/*" -type f -name "${pheno_line}_*BestMod.rds" 2>/dev/null)
      fi

      #if test model, get final model
      if [[ "$train_test" == "TEST" ]]; then
        mapfile -t matches < <(find "$save_dir" -path "*model_objs/*" -type f -name "${pheno_line}_*full_mod.rds" 2>/dev/null)
      fi

      #write
      if [ ${#matches[@]} -eq 1 ]; then
        mod="${matches[0]}"
        # tab-delimited: csv, model, traintest, split, total, save_path
        echo -e "$file\t$mod\t$train_test\t$total\t$split\t$save_path" >> "$config_file"
      elif [ ${#matches[@]} -eq 0 ]; then
        echo "Warning: No matching model found in '$save_dir' for '$pheno_line'" >&2
      else
        echo "Warning: Multiple matching files found in '$save_dir':" >&2
        printf '%s\n' "${matches[@]}" >&2
      fi
    done < "$pheno_list"

  done
  #add numbering (col 1 = slurm array task id)
  nl "$config_file" > temp.txt && mv temp.txt "$config_file"

  n=$(wc -l < "$config_file")
  echo "wrote $n rows to $config_file"
done
