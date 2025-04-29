#!/bin/bash
#script to write config file that will be used to models on 1/2 dataframe
#run from outside code dir

#PATHS
data_path=./data
config_path=./code/config_files
pheno_lists=./pheno_lists
log_scale="FALSE" #trying w/o scaling for now

if [ "$#" -eq 0 ]; then
  echo "No arguments provided. Need weighted TRUE/FALSE, total TRUE/FALSE"
  exit 1
  elif [ "$#" -ne 2 ]; then
  echo "Need weighted TRUE/FALSE, total TRUE/FALSE"
  exit 1
fi

echo "weighted = $1"
echo "include total value = $2"

#LOOP THROUGH 1/2 CSVS
for file in $(find $(realpath $data_path)  -type f -name "v3_CN_cleaned.csv") #for testing full df v3_CN_cleaned
do
  
  echo "prepping: $file"
  
  #get filename
  filename=$(basename -- "$file")
  filename="${filename%.*}"
  
  #LOOP THROUGH PHENO CATEGORIES
  for pheno_list in $(find $(realpath $pheno_lists) -type f -name "*.txt")
  do
    echo "pheno list: $pheno_list"
    
    #get filename -> freesurfer variable name
    pheno_cat=$(basename -- "$pheno_list")
    pheno_cat="${pheno_cat%.*}"
  
    #CREATE OUTPUT DIRS
    #make config file dir or remove old file if necessary
    config_file=$config_path/${filename}_${pheno_cat}_weight${1}_total${2}_scale${log_scale}config.txt
    if ! [ -d $config_path ]
    then
      mkdir $config_path
    elif [ -f $config_file ]
    then
      rm -rf $config_file
    fi
  
    touch $config_file
    
    #make output dir
      save_dir=./${filename}_train
      if ! [ -d $save_dir ]
      then
        mkdir $save_dir
      fi
      #name subdir based on whether total is controlled for
      if [[ $2 == "TRUE" ]]; then
        save_path=$save_dir/${pheno_cat}_total_mods
      elif [[ $2 == "FALSE" ]]; then
        save_path=$save_dir/${pheno_cat}_mods
      fi
      save_path=$(realpath $save_path) #get full paths
      if ! [ -d $save_path ]
      then
        mkdir $save_path
        mkdir $save_path/model_objs
        mkdir $save_path/centile_plots
        mkdir $save_path/worm_plots
        mkdir $save_path/cent_csvs
        mkdir $save_path/model_sums
      fi
      
      if [[ $pheno_cat == *"global"* ]]; then
        fs="fs_version_GM"
        tot="eTIV"
      elif [[ $pheno_cat == *"vols"* ]]; then
        fs="fs_version_GM"
        tot="TBV"
      elif [[ $pheno_cat == *"thickness"* ]]; then
        fs="fs_version_CT"
        tot="total.CT"
      elif [[ $pheno_cat == *"surf"* ]]; then
        fs="fs_version_SA"
        tot="total.SA"
      else
        echo "can't find appropriate variables"
      fi
    
      #LOOP THROUGH LAMBDAS
      for lambda in NULL #$(seq 100 100 50000)
      do
      
      if [[ $2 == "TRUE" ]]; then
      
        #LOOP THROUGH PHENOS
        while read -r pheno_line
        do
        # Write the CSV file path and the formula to the output file (tab-delimited)
          echo -e "$file\t$pheno_line\t$lambda\t$fs\t$tot\t$save_path\t$1\t$log_scale" >> "$config_file"
        done < "$pheno_list"
      
      elif [[ $2 == "FALSE" ]]; then
      
        #LOOP THROUGH PHENOS
        while read -r pheno_line
        do
        # Write the CSV file path and the formula to the output file (tab-delimited)
          echo -e "$file\t$pheno_line\t$lambda\t$fs\t$save_path\t$1\t$log_scale" >> "$config_file"
        done < "$pheno_list"
      fi
      
      done
      
   #add numbering
  nl "$config_file" > temp.txt && mv temp.txt "$config_file"
  done
done
