#!/bin/bash
#script to write config file that will be used to submit ComBatLS jobs using `subjobs_combatLS.sh`

base_path=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts

data_path=$base_path/data/to_combat #path to csv to be harmonized OR directory containing csvs
config_path=$base_path/code/config_files
pheno_path=$base_path/pheno_lists
save_path=$base_path/data/harmonized #path to save output csvs

#batch column name
batch="study_site"

#list covariate effects to preserve - MUST BE NUMERIC, CAN DUMMY-CODE
covar_list="logAge_days,sexMale,sexMale_x_logAge"

#combat configurations (formula for the above covars - keeping mu and sigma formulas identical for now)
mu_model="~pb(logAge_days)+sexMale+pb(sexMale_x_logAge)"
sigma_model="~pb(logAge_days)+sexMale+pb(sexMale_x_logAge)"

##########################################
####### DO NOT EDIT BELOW #######
##########################################
#make save path
if ! [ -d $save_path ]
	then
	mkdir $save_path
	mkdir $save_path/combat_objs
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
  
  # if data_path is csv, write out
	if [ -f $data_path ]
	then
		echo -e "$count\t$data_path\t$pheno_list\t$batch\t$covar_list\t$mu_model\t$sigma_model\t$full_save" >> "$config_path/combat_config.txt"
		count=$(( count+1 ))

	#else loop over data_path dir and match pheno name to csv name
	elif [ -d $data_path ]
	then
		#find corresponding data csvs
		for file in $(find $(realpath $data_path)  -type f -name $list_name*)
		do
		
		# Write the CSV file path and the formula to the output file (tab-delimited)
		echo -e "$count\t$file\t$pheno_list\t$batch\t$covar_list\t$mu_model\t$sigma_model\t$full_save" >> "$config_path/combat_config.txt"

		count=$(( count+1 ))
		
		done
	fi

done
