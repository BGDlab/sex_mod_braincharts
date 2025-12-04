#!/bin/bash
#script to write config file that will be used to test models on other sample
#run from outside code dir

#PATHS
data_path=./data
config_path=./code/config_files
pheno_lists=./pheno_lists

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --total) total="$2"; shift 2 ;;
    --log_age) log_age="$2"; shift 2 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

# Check that all required arguments are provided
if [[ -z "$total" || -z "$log_age" ]]; then
  echo "Missing arguments. Usage:"
  echo "--total TRUE/FALSE --log_age TRUE/FALSE"
  exit 1
fi

#LOOP THROUGH 1/2 CSVS
for split in A B 
do
  
  echo "prepping: $split"
  
    #make config file dir or remove old file if necessary
    config_file=$config_path/cv_sample_${split}_total${total}_logAge${log_age}_weighted_test_config.txt
    if ! [ -d $config_path ]
    then
      mkdir $config_path
    elif [ -f $config_file ]
    then
      rm -rf $config_file
    fi
  
    touch $config_file
  
  #get split df
  file=$(realpath $data_path/cv_sample_${split}_euler.csv)
  
  #LOOP THROUGH PHENO CATEGORIES
  for pheno_list in $(find $(realpath $pheno_lists) -type f -name "*.txt")
  do
    echo "pheno list: $pheno_list"

    #get filename -> total variable name
    pheno_cat=$(basename -- "$pheno_list")
    pheno_cat="${pheno_cat%.*}"
    
    if [[ "$total" == "FALSE" ]]; then
        tot="FALSE"
      elif [[ $pheno_cat == *"global"* ]]; then
        tot="TBV"
      elif [[ $pheno_cat == *"vols"* ]]; then
        tot="TBV"
      elif [[ $pheno_cat == *"thickness"* ]]; then
        tot="mean.CT"
      elif [[ $pheno_cat == *"surf"* ]]; then
        tot="total.SA"
      else
        echo "can't find appropriate variables"
      fi
  
    #CREATE OUTPUT DIRS
    #make output dir
      save_dir=./cv_sample_${split}_test
      if ! [ -d $save_dir ]
      then
        mkdir $save_dir
      fi
      save_path=$save_dir/weighted_${pheno_cat}_total${total}_logAge${log_age}_pbmods
      #save_path=$(realpath $save_path) #get full paths
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
      
      #SET SEARCH PATH FOR TRAINING MODELS
      #replace 'test' with 'train' and swap A and B
      search_path=$(echo "$save_dir/weighted_${pheno_cat}_total${total}_*logAge${log_age}_pbmods" | \
        sed -E 's/test/train/g; s/_A_/_TEMP_/g; s/_B_/_A_/g; s/_TEMP_/_B_/g') #wildcard to allow for old logPheno naming convention
      
      #LOOP THROUGH BESTMODS
      while read -r pheno_line
      do
        
        # find training BestModel from model selection
        mapfile -t matches < <(find $search_path -type f -name "${pheno_line}_*weighted.rds")
        if [ ${#matches[@]} -eq 1 ]; then
          og_mod=$(realpath "${matches[0]}")
          
          # Write the CSV file path and the formula to the output file (tab-delimited)
          echo -e "$file\t$og_mod\t$save_path\t$tot" >> "$config_file"
        elif [ ${#matches[@]} -eq 0 ]; then
          echo "Warning: No matching file found in '$search_path' for prefix '$pheno_line' and suffix 'weighted.rds'" >&2
          #exit 1
        else
          echo "Warning: Multiple matching files found in '$search_path':" >&2
          printf '%s\n' "${matches[@]}" >&2
          #exit 1
        fi
      done < "$pheno_list"

  done
   #add numbering
  nl "$config_file" > temp.txt && mv temp.txt "$config_file"
done
