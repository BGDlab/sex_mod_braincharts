#!/bin/bash
#script to write config file that will be used to submit ComBatLS jobs using `subjobs_combatLS.sh`
#run from outside code dir

#list csvs to be harmonized
data_path=./data/to_combat
config_path=./code/config_files
pheno_path=./pheno_lists
save_path=./data/harmonized

#batch
batch="study_site"

#list covariate effects to preserve
covar_list="logAge_days,sexMale,sexMale_x_logAge"

#combat configurations (formula for the above covars - keeping mu and sigma formulas identical for now)
#knots from df4 - will need to update manually to change
combat_mod="~ns(logAge_days,knots=c(2.7645,3.3654,3.966))+sexMale+ns(sexMale_x_logAge,knots=c(2.765,3.3654,3.966))"

#make save path
if ! [ -d $save_path ]
	then
	mkdir $save_path
fi
full_save=$(realpath $save_path)

#make config file dir or remove old file if necessary
if ! [ -d $config_path ]
	then
	mkdir $config_path
	elif [ -f $config_path/combat_config.txt ]
	then
	rm -rf $config_path/combat_config.txt
fi

# create the output file
touch $config_path/combat_config.txt

count=1

#loop through each list of phenos
for pheno_list in $(find $(realpath $pheno_path)  -type f -name "*.rds")
do

  list_name=$(basename $pheno_list .rds)
  echo "pheno_list: $list_name"

  #find corresponding data csvs
  for file in $(find $(realpath $data_path)  -type f -name $list_name*)
  do
  
  # Write the CSV file path and the formula to the output file (tab-delimited)
  echo -e "$count\t$file\t$pheno_list\t$batch\t$covar_list\t$combat_mod\t$combat_mod\t$full_save" >> "$config_path/combat_config.txt"

  count=$(( count+1 ))
  
  done
done
