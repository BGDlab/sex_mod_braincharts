#!/bin/bash
#script to write config file that will be used to submit ComBatLS jobs using `subjobs_combatLS.sh`
#run from outside code dir

#list csvs to be harmonized
data_path=./data/to_combat
config_path=./code/config_files

#batch
batch="study"

#list covariate effects to preserve
covar_list="logAge_days,sexMale,sexMale_x_logAge"

#list desired combat configurations (formula for the above covars - keeping mu and sigma formulas identical for now)
combat_list='"~pb(logAge_days, method='GAIC', k=log(nrow(df)) + sexMale + pb(sexMale_x_logAge, method='GAIC', k=log(nrow(df))" \
"~ ns(logAge_days, df=30) + sexMale + ns(sexMale_x_logAge, df=10)" \
"~ ns(logAge_days, df=20) + sexMale + ns(sexMale_x_logAge, df=5)" \
"~ ns(logAge_days, knots=c(2.562590, 3.039711, 3.407688, 3.766710, 3.960530, 4.215803, 4.407688)) + sexMale + ns(sexMale_x_logAge, df=10)"'

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


#loop through each CSV of IDPs separated by phenotype category
for file in $(find $(realpath $data_path)  -type f )
do
  echo "prepping: $file"
  for combat_config in $combat_list
  do
  # Write the CSV file path and the formula to the output file (tab-delimited)
  echo -e "$file\t$covar_list\t$combat_config\t$combat_config" >> "$config_path/combat_config.txt"
  done
done
