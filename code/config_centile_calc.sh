#!/bin/bash
#script to write config file that will be used to test models on other sample
#run from outside code dir

#PATHS
data_path=./data
config_path=./code/config_files
pheno_lists=./pheno_lists

rerun="FALSE"
age2plus="FALSE"
weighted="FALSE"

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --total) total="$2"; shift 2 ;;
    --log_age) log_age="$2"; shift 2 ;;
    --rerun) rerun="$2"; shift 2 ;;
    --age2plus) age2plus="$2"; shift 2 ;;
    --weighted) weighted="$2"; shift 2 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

# Check that all required arguments are provided
if [[ -z "$total" || -z "$log_age" ]]; then
  echo "Missing arguments. Usage:"
  echo "--total TRUE/FALSE --log_age TRUE/FALSE [--rerun TRUE/FALSE] [--age2plus TRUE/FALSE] [--weighted TRUE/FALSE]"
  exit 1
fi

if [[ "$age2plus" == "TRUE" && "$weighted" == "TRUE" ]]; then
  echo "Error: --age2plus and --weighted cannot both be TRUE (no combined dirs exist)"
  exit 1
fi

# Print arguments
echo "include total value = $total"
echo "log scale age = $log_age"
echo "age2plus mode = $age2plus"
echo "weighted mode = $weighted"

#LOOP THROUGH 1/2 CSVS
for split in A B
do
  
  echo "prepping: $split"
  
    #make config file dir or remove old file if necessary
    if [[ "$age2plus" == "TRUE" ]]; then
      mode_tag="age2plus_centtest"
    elif [[ "$weighted" == "TRUE" ]]; then
      mode_tag="weighted_centtest"
    else
      mode_tag="centtest"
    fi
    if [[ "$rerun" == "TRUE" ]]; then
      date_tag=$(date +%Y%m%d)
      config_file=$config_path/cv_sample_${split}_total${total}_logAge${log_age}_${mode_tag}_rerun${date_tag}_config.txt
    else
      config_file=$config_path/cv_sample_${split}_total${total}_logAge${log_age}_${mode_tag}_config.txt
    fi
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
  
    #CREATE OUTPUT DIRS
    
    #make output dir
      save_dir=./cv_sample_${split}_test

      if [[ "$age2plus" == "TRUE" ]]; then
        save_path=$save_dir/age2plus_${pheno_cat}_total${total}_logAge${log_age}_pbmods
      elif [[ "$weighted" == "TRUE" ]]; then
        save_path=$save_dir/weighted_${pheno_cat}_total${total}_logAge${log_age}_pbmods
      else
        save_path=$save_dir/${pheno_cat}_total${total}_logAge${log_age}_pbmods
      fi
      save_path=$(realpath $save_path)

      #LOOP THROUGH TEST MODELS
      while read -r pheno_line
      do
      
      #LOOP THROUGH DX
      for dx in "SCZ" "ALZ" "ASD" "MDD" "GAD" "ADHD"
      do
      
      #skip already tested models if rerunning
      if [[ "$rerun" == "TRUE" ]]; then
          test_csv=$(ls ${save_path}/cent_csvs/${pheno_line}*PT_${dx}_cent.csv 2>/dev/null | head -n 1 || true)
          if [[ -n "$test_csv" ]]; then
            echo "Skipping $pheno_line (centiles found)"
            continue
          fi
      fi
      
        #find csv to test in - handle optional _logPheno*_ in filename
        mapfile -t pt_file_matches < <(find $(realpath $data_path/case_control_dfs) -type f -name "cv_sample_${split}_${dx}.csv" 2>/dev/null)
        if [ ${#pt_file_matches[@]} -gt 1 ]; then
          echo "Error: Multiple CSV files found for $pheno_line:"
          printf '%s\n' "${pt_file_matches[@]}"
          exit 1
        elif [ ${#pt_file_matches[@]} -eq 0 ]; then
          echo "Warning: No CSV found for $pheno_line, skipping"
          continue
        else
          pt_file="${pt_file_matches[0]}"
        fi
        
        #find csv model test was fit in
        if [[ "$age2plus" == "TRUE" ]]; then
          mapfile -t file_matches < <(find $(realpath $data_path/cv_sample_${split}_age2plus) -type f -name "${pheno_line}_total${total}_test_df.csv" 2>/dev/null)
        elif [[ "$weighted" == "TRUE" ]]; then
          mapfile -t file_matches < <(find $(realpath $data_path/cv_sample_${split}_euler) -type f -name "${pheno_line}_total${total}_test_df.csv" 2>/dev/null)
        else
          mapfile -t file_matches < <(find $(realpath $data_path/cv_sample_${split}_dfs) -type f -name "${pheno_line}_total${total}*logAge${log_age}.csv" 2>/dev/null)
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

	    #find test model - handle optional _logPheno*_ in directory names
        if [[ "$age2plus" == "TRUE" ]]; then
          mapfile -t matches < <(
            find "$(realpath "$save_dir")" \
              -path "*age2plus_${pheno_cat}_total${total}*logAge${log_age}_pbmods/model_objs/*" \
              -type f \
              -name "${pheno_line}_*full_mod.rds" \
              ! -path "*weighted*" \
              2>/dev/null
          )
        elif [[ "$weighted" == "TRUE" ]]; then
          mapfile -t matches < <(
            find "$(realpath "$save_dir")" \
              -path "*weighted_${pheno_cat}_total${total}*logAge${log_age}_pbmods/model_objs/*" \
              -type f \
              -name "${pheno_line}_*full_mod.rds" \
              ! -path "*age2plus*" \
              2>/dev/null
          )
        else
          mapfile -t matches < <(
            find "$(realpath "$save_dir")" \
              -path "*${pheno_cat}_total${total}*logAge${log_age}_pbmods/model_objs/*" \
              -type f \
              -name "${pheno_line}_*full_mod.rds" \
              ! -path "*weighted*" \
              ! -path "*age2plus*" \
              2>/dev/null
          )
        fi
        if [ ${#matches[@]} -eq 1 ]; then
          og_mod="${matches[0]}"
          
          # Write the CSV file path and the formula to the output file (tab-delimited)
          echo -e "$pt_file\t$file\t$og_mod\t$save_path\t$dx" >> "$config_file"
        elif [ ${#matches[@]} -eq 0 ]; then
          echo "Warning: No matching file found in '$save_dir' for prefix '$pheno_line' and suffix 'full_mod.rds'" >&2
          #exit 1
        else
          echo "Warning: Multiple matching files found in '$save_dir':" >&2
          printf '%s\n' "${matches[@]}" >&2
          #exit 1
        fi
        done
      done < "$pheno_list"

  done
   #add numbering
  nl "$config_file" > temp.txt && mv temp.txt "$config_file"
done
