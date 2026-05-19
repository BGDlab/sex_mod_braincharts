#!/bin/bash
#script to write config file that will be used to test models on other sample
#run from outside code dir

#PATHS
data_path=./data
config_path=./code/config_files
save_path=./dx_tests

rerun="FALSE"

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --total) total="$2"; shift 2 ;;
    --log_age) log_age="$2"; shift 2 ;;
    --rerun) rerun="$2"; shift 2 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

# Check that all required arguments are provided
if [[ -z "$total" || -z "$log_age" ]]; then
  echo "Missing arguments. Usage:"
  echo "--total TRUE/FALSE --log_age TRUE/FALSE [--rerun TRUE/FALSE]"
  exit 1
fi

# Print arguments
echo "include total value = $total"
echo "log scale age = $log_age"

#make config file dir or remove old file if necessary
if [[ "$rerun" == "TRUE" ]]; then
  date_tag=$(date +%Y%m%d)
  config_file=$config_path/cv_mods_total${total}_logAge${log_age}_centext_rerun${date_tag}_config.txt
  else
  config_file=$config_path/cv_mods_total${total}_logAge${log_age}_centext_config.txt
  fi
if ! [ -d $config_path ]
then
  mkdir $config_path
  elif [ -f $config_file ]
  then
    rm -rf $config_file
  fi
  
  touch $config_file

#make output dir
if ! [ -d $save_path ]
then
mkdir $save_path
fi
    
#LOOP THROUGH DX
for dx in "SCZ" "ALZ" "ASD" "MDD" "GAD" "ADHD"
do
      
  #skip already tested models if rerunning
  if [[ "$rerun" == "TRUE" ]]; then
    test_csv=$(ls ${save_path}/cent_csvs/${dx}_extcent_diffs.csv 2>/dev/null | head -n 1 || true)
    if [[ -n "$test_csv" ]]; then
      echo "Skipping $dx (extreme diffs found)"
    continue
      fi
    fi
      
  #write config file
  echo -e "$dx\t$total" >> "$config_file"

done
 #add numbering
  nl "$config_file" > temp.txt && mv temp.txt "$config_file"
