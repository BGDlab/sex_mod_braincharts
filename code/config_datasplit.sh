#!/bin/bash
#script to write config file that will be used to submit pre-ComBat data-splitting jobs using `subjobs_datasplit.sh`
#run from outside code dir

#list csvs to be harmonized
data_path=./data/csvs_by_pheno
config_path=./code/config_files
pheno_lists=./pheno_lists
save_path=./data/to_combat

#make config file dir or remove old file if necessary
if ! [ -d $config_path ]
	then
	mkdir $config_path
	elif [ -f $config_path/datasplit_config.txt ]
	then
	rm -rf $config_path/datasplit_config.txt
fi

#make output dir
if ! [ -d $save_path ]
	then
	mkdir $save_path
fi
save_path=$(realpath $save_path)

# create the output file
touch $config_path/datasplit_config.txt

count=1

#loop through each CSV of IDPs separated by phenotype category
for file in $(find $(realpath $data_path)  -type f -name "*.csv")
do
  echo "prepping: $file"
  #find corresponding pheno list
  filename=$(basename -- "$file")
  filename="${filename%.*}"
  pheno_list=$(find $(realpath $pheno_lists) -type f -name "$filename.rds")
  echo "pheno list: $pheno_list"

  # Write the CSV file path and the formula to the output file (tab-delimited)
  echo -e "$count\t$file\t$pheno_list\t$save_path\t$filename" >> "$config_path/datasplit_config.txt"
  
  count=$(( count+1 ))
done
