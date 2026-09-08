#!/bin/bash

#run from outside code dir

#PATHS
data_path=./data
config_path=./code/config_files
pheno_lists=./pheno_lists

#make config file dir or remove old file if necessary
config_file=$config_path/cv_sample_logAgeTRUE_clean_config.txt

if ! [ -d $config_path ]
then
  mkdir $config_path
elif [ -f $config_file ]
then
  rm -rf $config_file
fi

touch $config_file

#make save dir
save_dir=./models_to_share
if ! [ -d $save_dir ]
then
  mkdir $save_dir
fi
save_path=$(realpath $save_dir)

#fig outputs
fig_save_dir=./sharing_code/figs
if ! [ -d $fig_save_dir ]
then
  mkdir $fig_save_dir
fi
fig_save_dir=$(realpath $fig_save_dir)

#LOOP THROUGH SPLITS
for split in A B
do

  echo "prepping: $split"
  
  #LOOP THROUGH TOTAL TRUE/FALSE
  for total in TRUE FALSE
  do
  echo "total-size-corrected = $total"

  #LOOP THROUGH PHENO CATEGORIES
  for pheno_list in $(find $(realpath $pheno_lists) -type f -name "*.txt")
  do
    echo "pheno list: $pheno_list"

    pheno_cat=$(basename -- "$pheno_list")

    #LOOP THROUGH PHENOS
    while read -r pheno_line
    do

      #find csv - handle optional _logPheno*_ in filename
      #if pheno_cat == cortical_thickness & total == TRUE -> age2plus
      if [[ "$pheno_cat" == "cortical_thickness" && "$total" == "TRUE" ]]; then
        model_path=./cv_sample_${split}_test/age2plus_${pheno_cat}_total${total}_logAgeTRUE_pbmods
        mapfile -t file_matches < <(find $(realpath $data_path/cv_sample_${split}_age2plus) -type f -name "${pheno_line}_total${total}*logAgeTRUE.csv" 2>/dev/null)
      else
        model_path=./cv_sample_${split}_test/${pheno_cat}_total${total}_logAgeTRUE_pbmods
        mapfile -t file_matches < <(find $(realpath $data_path/cv_sample_${split}_dfs) -type f -name "${pheno_line}_total${total}*logAgeTRUE.csv" 2>/dev/null)
      fi
      
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
      
      #get model
      mapfile -t matches < <(find "$model_path" -path "*model_objs/*" -type f -name "${pheno_line}_*full_mod.rds" 2>/dev/null)

      #write
      if [ ${#matches[@]} -eq 1 ]; then
        mod="${matches[0]}"
        # tab-delimited: csv, model, traintest, split, total, save_path
        echo -e "$file\t$mod\t$pheno_cat\t$total\t$split\t$save_path\t$fig_save_dir" >> "$config_file"
      elif [ ${#matches[@]} -eq 0 ]; then
        echo "Warning: No matching model found in '$model_path' for '$pheno_line'" >&2
      else
        echo "Warning: Multiple matching files found in '$model_path':" >&2
        printf '%s\n' "${matches[@]}" >&2
      fi
    done < "$pheno_list"
  done
  done
done

#add numbering (col 1 = slurm array task id)
nl "$config_file" > temp.txt && mv temp.txt "$config_file"

n=$(wc -l < "$config_file")
echo "wrote $n rows to $config_file"

