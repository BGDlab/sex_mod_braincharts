#!/bin/bash
#script to write config file for bootstrapping test models (one line per split x pheno)
#each line is submitted as its own array job by submit_boot_mods.sh (one task per bootstrap)
#run from outside code dir

#PATHS
data_path=./data
config_path=./code/config_files
#hardcoding sim df for consistency across phenos and splits
sim_df=$(realpath $data_path/v3_full_sim.rds)

pheno_list=./pheno_lists/global_vols.txt
seed=20261021
strata="study_site"

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --total) total="$2"; shift 2 ;;
    --log_age) log_age="$2"; shift 2 ;;
    --pheno_list) pheno_list="$2"; shift 2 ;;
    --seed) seed="$2"; shift 2 ;;
    --strata) strata="$2"; shift 2 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

# Check that all required arguments are provided
if [[ -z "$total" || -z "$log_age" ]]; then
  echo "Missing arguments. Usage:"
  echo "--total TRUE/FALSE --log_age TRUE/FALSE [--pheno_list path] [--seed int] [--strata study/study_site]"
  exit 1
fi

echo "include total value = $total"
echo "log scale age = $log_age"
echo "pheno list = $pheno_list"
echo "seed = $seed"
echo "strata = $strata x sex"

pheno_list=$(realpath $pheno_list)
pheno_cat=$(basename -- "$pheno_list")
pheno_cat="${pheno_cat%.*}"

#make config file dir or remove old file if necessary
config_file=$config_path/boot_${pheno_cat}_total${total}_logAge${log_age}_config.txt
if ! [ -d $config_path ]
then
  mkdir $config_path
elif [ -f $config_file ]
then
  rm -rf $config_file
fi
touch $config_file

#LOOP THROUGH SPLITS
for split in A B
do
  echo "prepping: $split"

  test_dir=./cv_sample_${split}_test
  #training models come from the other split
  train_dir=$(echo "$test_dir" | \
    sed -E 's/test/train/g; s/_A_/_TEMP_/g; s/_B_/_A_/g; s/_TEMP_/_B_/g')

  while read -r pheno_line
  do
    #TBV is not scored as a total-corrected pheno
    if [[ "$total" == "TRUE" && "$pheno_line" == "TBV" ]]; then
      echo "Skipping TBV for total=TRUE"
      continue
    fi

    #csv the test model was fit on - handle optional _logPheno*_ in filename
    mapfile -t file_matches < <(find $(realpath $data_path/cv_sample_${split}_dfs) -type f -name "${pheno_line}_total${total}*logAge${log_age}.csv" 2>/dev/null)
    if [ ${#file_matches[@]} -ne 1 ]; then
      echo "Warning: ${#file_matches[@]} CSVs found for $pheno_line in split $split, skipping" >&2
      printf '%s\n' "${file_matches[@]}" >&2
      continue
    fi
    file="${file_matches[0]}"

    #training BestMod (fixed lambdas)
    mapfile -t train_matches < <(find "$(realpath "$train_dir")" -path "*${pheno_cat}_total${total}*logAge${log_age}_pbmods/model_objs/*" -type f -name "${pheno_line}_*BestMod.rds" 2>/dev/null)
    if [ ${#train_matches[@]} -ne 1 ]; then
      echo "Warning: ${#train_matches[@]} training BestMods found for $pheno_line in $train_dir, skipping" >&2
      printf '%s\n' "${train_matches[@]}" >&2
      continue
    fi
    train_mod="${train_matches[0]}"

    #test full model (warm start) - exclude weighted & age2plus versions
    mapfile -t test_matches < <(find "$(realpath "$test_dir")" -path "*${pheno_cat}_total${total}*logAge${log_age}_pbmods/model_objs/*" \
      -not -path "*weight*" -not -path "*age2plus*" -type f -name "${pheno_line}_*test_full_mod.rds" 2>/dev/null)
    if [ ${#test_matches[@]} -ne 1 ]; then
      echo "Warning: ${#test_matches[@]} test full_mods found for $pheno_line in $test_dir, skipping" >&2
      printf '%s\n' "${test_matches[@]}" >&2
      continue
    fi
    test_mod="${test_matches[0]}"

    #make output dir
    save_path=./cv_sample_${split}_boot/${pheno_cat}_total${total}_logAge${log_age}_boot/${pheno_line}
    mkdir -p $save_path
    save_path=$(realpath $save_path)

    echo -e "$file\t$train_mod\t$test_mod\t$save_path\t$sim_df\t$seed\t$strata" >> "$config_file"
  done < "$pheno_list"
done

#add numbering
nl "$config_file" > temp.txt && mv temp.txt "$config_file"
echo "wrote $(wc -l < $config_file) lines to $config_file"

#open permissions so colleagues can submit/rerun jobs and write outputs & logs here
#(chmod only works on files you own - skip anything a colleague created)
mkdir -p ./code/jobfiles/boot ./code/jobfiles/combine_boot
for split in A B
do
  [ -d ./cv_sample_${split}_boot ] && chmod -R 777 ./cv_sample_${split}_boot 2>/dev/null
done
chmod 777 $config_path $config_file ./code/jobfiles/boot ./code/jobfiles/combine_boot 2>/dev/null
true
