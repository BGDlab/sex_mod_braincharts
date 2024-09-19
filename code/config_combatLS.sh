#script to write config file that will be used to submit ComBatLS jobs using `subjobs_combatLS.sh`

#list csvs to be harmonized
data_path=../data
config_path=./config_files

#batch
batch="study"

#list covariate effects to preserve
covar_list=

#list desired combat configurations (formula for the above covars)
combat_list=("y=mx+b" "y ~ pb(x) + b")

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
for file in $(find "$data_path" -type f \( -name "*gmv.csv" -o -name "*sa.csv" -o -name "*ct.csv" -o -name "*volglob.csv" -o -name "*sub.csv"\))
do
  echo "prepping: $file"
  for combat_config in combat_list
  do
  # Write the CSV file path and the formula to the output file (tab-delimited)
  echo -e "$file\t$combat_config" >> "$config_path/combat_config.txt"
  done
done