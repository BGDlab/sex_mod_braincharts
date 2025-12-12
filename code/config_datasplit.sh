#!/bin/bash
#script to write config file that will be used to submit pre-ComBat data-splitting jobs using `subjobs_datasplit.sh`
#run from outside code dir

# Set values for pheno_lists
base_path=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts
pheno_lists=$base_path/pheno_lists
config_path=$base_path/code/config_files
file=$base_path/data/v3_CN_cleaned.csv
save_path=$base_path/data/to_combat
batch_arg="study_site"

file=$(realpath $file)

#make config file dir or remove old file if necessary
if ! [ -d $config_path ]
	then
	mkdir $config_path
	elif [ -f $config_path/datasplit_config.txt ]
	then
	rm -rf $config_path/datasplit_config.txt
fi

# create the output file
tmp_path=$(realpath $config_path)
touch $config_path/datasplit_config.txt

#make output dir
if ! [ -d $save_path ]
	then
	mkdir $save_path
fi

count=1

#loop through each list of IDPs separated by phenotype category
for pheno_list in $(find $(realpath $pheno_lists)  -type f -name "*.rds")
do
  echo "prepping: $pheno_list"

  # write out args
  echo -e "$count\t$file\t$pheno_list\t$save_path\t$batch_arg" >> "$config_path/datasplit_config.txt"
  
  count=$(( count+1 ))
done
