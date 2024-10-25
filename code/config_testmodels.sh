#!/bin/bash
#script to write config file that will be used to submit test models with varying smooths
#run from outside code dir

#paths
data_path=./data/csvs_by_pheno
config_path=./code/config_files
pheno_lists=./pheno_lists
knot_lists=./code/knot_lists.RDS

######################################
#Configuration Options
######################################
if [ $1 = "no_fs" ]
then
  config_file=$config_path/testmodels_config_nofs.txt
  save_path=./test_sex_mod_nofs
  log_scale="FALSE"
  fs="FALSE"
  elif [ $1 = "log" ]
  then
  config_file=$config_path/testmodels_config_log.txt
  save_path=./test_sex_mod_log
  log_scale="TRUE"
  fs="TRUE"
  elif [ $1 = "basic" ]
  then
  config_file=$config_path/testmodels_config.txt
  save_path=./test_sex_mod
  log_scale="FALSE"
  fs="TRUE"
fi
######################################

#get full paths
save_path=$(realpath $save_path)
knot_lists=$(realpath $knot_lists)

#make config file dir or remove old file if necessary
if ! [ -d $config_path ]
	then
	mkdir $config_path
	elif [ -f $config_file ]
	then
	rm -rf $config_file
fi

#make output dir
if ! [ -d $save_path ]
	then
	mkdir $save_path
	mkdir $save_path/model_objs
	mkdir $save_path/centile_plots
	mkdir $save_path/worm_plots
	mkdir $save_path/cent_csvs
fi

# create the output file
touch $config_file

#loop through each CSV of IDPs separated by phenotype category
for file in $(find $(realpath $data_path)  -type f -name "*.csv")
do
  echo "prepping: $file"
  #find corresponding pheno list
  filename=$(basename -- "$file")
  filename="${filename%.*}"
  pheno_list=$(find $(realpath $pheno_lists) -type f -name "$filename.txt")
  echo "pheno list: $pheno_list"
  
  while read -r pheno_line
  do

  # Write the CSV file path and the formula to the output file (tab-delimited)
  echo -e "$file\t$pheno_line\t$knot_lists\t$save_path\t$log_scale\t$fs" >> "$config_file"
  done < "$pheno_list"
done

#add numbering
nl "$config_file" > temp.txt && mv temp.txt "$config_file"
