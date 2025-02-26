#!/bin/bash
#script to write config file that will be used to models on 1/2 dataframe
#run from outside code dir

#CHECK ARG
if [[ $# -ne 1 ]]; then
    echo 'pass arg for moments to include fs covariate: "both", "mu", or "none"' >&2
    exit 1
fi

case $1 in
    both|mu|none)  # Ok
        ;;
    *)
        # The wrong first argument.
        echo 'pass arg for moments to include fs covariate: "both", "mu", or "none"' >&2
        exit 1
esac

#PATHS
data_path=./data
config_path=./code/config_files
pheno_lists=./pheno_lists
log_scale="FALSE" #trying w/o scaling for now

#LOOP THROUGH 1/2 CSVS
for file in $(find $(realpath $data_path)  -type f -name "cv_sample*.csv")
do
  
  echo "prepping: $file"
  
  #get filename
  filename=$(basename -- "$file")
  filename="${filename%.*}"
  
  #CREATE OUTPUT DIRS
  touch $config_file
  #make config file dir or remove old file if necessary
  config_file=$config_path/$filename_$1_config.txt
  if ! [ -d $config_path ]
  then
    mkdir $config_path
  elif [ -f $config_file ]
  then
    rm -rf $config_file
  fi
  
  #make output dir
    save_path=./$filename_$1_mods
    save_path=$(realpath $save_path) #get full paths
    if ! [ -d $save_path ]
    then
      mkdir $save_path
      mkdir $save_path/model_objs
      mkdir $save_path/centile_plots
      mkdir $save_path/worm_plots
      mkdir $save_path/cent_csvs
      mkdir $save_path/model_sums #NEED TO IMPLEMENT
      mkdir $save_path/z_csvs #NEED TO IMPLEMENT
    fi
  
  #LOOP THROUGH PHENO CATEGORIES
  for pheno_list in $(find $(realpath $data_path) -type f -name "*.txt")
  do
    echo "pheno list: $pheno_list"
    
    #get filename -> freesurfer variable name
    pheno_cat=$(basename -- "$pheno_list")
    pheno_cat="${pheno_cat%.*}"
    
    if [[ $pheno_cat == *"vols"* ]]; then
      fs="freesurfer_version_GM"
    elif [[ $pheno_cat == *"thickness"* ]]; then
      fs="freesurfer_version_CT"
    elif [[ $pheno_cat == *"surf"* ]]; then
      fs="freesurfer_version_SA"
    else
      echo "can't find appropriate fs version"
    fi
    
    #LOOP THROUGH PHENOS
    while read -r pheno_line
    do
    # Write the CSV file path and the formula to the output file (tab-delimited)
      echo -e "$file\t$pheno_line\t$fs\t$1\t$save_path\t$log_scale" >> "$config_file"
    done < "$pheno_list"
  done
  
  #add numbering
  nl "$config_file" > temp.txt && mv temp.txt "$config_file"
done
