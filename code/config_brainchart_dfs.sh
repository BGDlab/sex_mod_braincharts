#!/bin/bash
#script to write config file that will prep dataframe for each pheno
#run from outside code dir!!!!

#PATHS
data_path=./data
config_path=./code/config_files
pheno_lists=./pheno_lists
file=./data/v3_CN_cleaned.csv

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --log_pheno) log_pheno="$2"; shift 2 ;;
    --total) total="$2"; shift 2 ;;
    --log_age) log_age="$2"; shift 2 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

# Check that all required arguments are provided
if [[ -z "$log_pheno" || -z "$total" || -z "$log_age" ]]; then
  echo "Missing arguments. Usage:"
  echo "--log_pheno TRUE/FALSE --total TRUE/FALSE --log_age TRUE/FALSE"
  exit 1
fi

# Print arguments
echo "log scale pheno = $log_pheno"
echo "include total value = $total"
echo "log scale age = $log_age"

  #make new subdir
  save_path=$data_path/pheno_dfs_total${total}
  if ! [ -d $save_path ]
    then
    mkdir $save_path
  fi
  save_path=$(realpath $save_path)
  
  #CREATE OUTPUT DIRS
  #make config file dir or remove old file if necessary
  config_file=$config_path/braincharts_logPheno${log_pheno}_total${total}_logAge${log_age}_df_config.txt
  if ! [ -d $config_path ]
  then
    mkdir $config_path
  elif [ -f $config_file ]
  then
    rm -rf $config_file
  fi
  
  #LOOP THROUGH PHENO CATEGORIES
  for pheno_list in $(find $(realpath $pheno_lists) -type f -name "*.txt")
  do
    echo "pheno list: $pheno_list"
    
    #get filename -> freesurfer variable name
    pheno_cat=$(basename -- "$pheno_list")
    pheno_cat="${pheno_cat%.*}"
  
    touch $config_file
      
    if [[ $pheno_cat == *"global"* ]]; then
        fs="fs_version_GM"
        tot="TBV"
    elif [[ $pheno_cat == *"vols"* ]]; then
        fs="fs_version_GM"
        tot="TBV"
    elif [[ $pheno_cat == *"thickness"* ]]; then
        fs="fs_version_CT"
        tot="mean.CT"
    elif [[ $pheno_cat == *"surf"* ]]; then
        fs="fs_version_SA"
        tot="total.SA"
    else
        echo "can't find appropriate variables"
    fi
    
    if [[ $total == "TRUE" ]]; then
      
      #LOOP THROUGH PHENOS
      while read -r pheno_line
      do
      # Write the CSV file path and the formula to the output file (tab-delimited) test
        echo -e "$file\t$pheno_line\t$fs\t$tot\t$log_pheno\t$log_age\t$save_path" >> "$config_file"
      done < "$pheno_list"
      
    elif [[ $total == "FALSE" ]]; then
      
      #LOOP THROUGH PHENOS
      while read -r pheno_line
      do
      # Write the CSV file path and the formula to the output file (tab-delimited)
        echo -e "$file\t$pheno_line\t$fs\tNULL\t$log_pheno\t$log_age\t$save_path" >> "$config_file"
      done < "$pheno_list"
    fi
  done
  #add numbering
  nl "$config_file" > temp.txt && mv temp.txt "$config_file"


chmod 777 -R $save_path
