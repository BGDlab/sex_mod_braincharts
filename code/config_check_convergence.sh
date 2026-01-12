#!/bin/bash
#run from outside code dir!!!!

#PATHS
dir_path=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts
config_path=./code/config_files
pheno_lists=./pheno_lists


config_file=$config_path/converge_check_config.txt
  if ! [ -d $config_path ]
  then
    mkdir $config_path
  elif [ -f $config_file ]
  then
    rm -rf $config_file
  fi
  
for dir in $(find "$dir_path" -type d -path "*/cv_sample_?_train")
do
  #LOOP THROUGH PHENO CATEGORIES
  for pheno_list in $(find $(realpath $pheno_lists) -type f -name "*.txt")
  do
  #LOOP THROUGH PHENOS
      while read -r pheno_line
      do
      echo -e "$dir\t$pheno_line" >> "$config_file"
      done < "$pheno_list"

  done
done
#add numbering
nl "$config_file" > temp.txt && mv temp.txt "$config_file"
