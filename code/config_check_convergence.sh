#!/bin/bash
#run from outside code dir!!!!

#PATHS
dir_path=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts
config_path=./code/config_files

config_file=$config_path/converge_check_config.txt
  if ! [ -d $config_path ]
  then
    mkdir $config_path
  elif [ -f $config_file ]
  then
    rm -rf $config_file
  fi
  
  find "$dir_path" -type d -path "*/cv_sample_?_train/*/model_objs" > "$config_file"

  #add numbering
  nl "$config_file" > temp.txt && mv temp.txt "$config_file"
