#!/bin/bash
# Apply ComFam() with different arguments to specified dataset

#STEP 2 OF ANALYSIS PIPELINE#

#######################################################################################
# SET PATHS
img=/cbica/home/gardnerm/software/containers/r_gamlss_0.0.1.sif #singularity image
base=/cbica/home/gardnerm/combat-biovar #base path (cubic)
save_path=$base/data/ukb_to_model
#######################################################################################
helpFunction()
{
   echo ""
   echo "Usage: $0 -c data.csv -t true -b 'batch variable'"
   echo -e "\t-c Path to data as .csv file"
   echo -e "\t-t true/false to log-transform"
   echo -e "\t-t name of variable (e.g. site, study) that defines batches"
   exit 1 # Exit script after printing help
}
#######################################################################################
# GET ARGS
while getopts ":c:t:b:" opt
do
   case "$opt" in
	c ) csv="$OPTARG" ;;
	t ) transform="$OPTARG" ;;
	b ) batch="$OPTARG" ;;
   esac
done

if [ -z "$csv" ] || [ -z "$transform" ] || [ -z "$batch" ]
then
   echo "provide: path to data and whether to log-transform global vols before running ComBat";
   helpFunction
fi
#######################################################################################
# DEFINE COMBAT .R SCRIPT (DEPENDING ON TRANSFORMATION OPTION)
if [ "${transform,,}" = true ]
then
cf_script=$base/cubic_scripts/R_scripts/combat_apply_w_transform.R #path to .R script - set to apply log-transform to global vols before combatting
elif [ "${transform,,}" = false ]
then
cf_script=$base/cubic_scripts/R_scripts/combat_apply.R #path to .R script
else
echo "Help! Invalid selection for transform true/false"
fi
#######################################################################################
cd $base/cubic_scripts #to source functions correctly
#######################################################################################
# MAKE DIRECTORIES 
#study dir
study_dir=$base/ukb_basic
if ! [ -d $study_dir ]
	then
	mkdir $study_dir
	fi

#qsub script & outputs
bash_dir=$study_dir/combat_qsub_scripts
if ! [ -d $bash_dir ]
	then
	mkdir $bash_dir
	#else
	#remove old error messages if necessary
	#rm -rf $bash_dir/*.txt
	fi

#combat output obj dir
if ! [ -d ${save_path}/combat_objs ]
	then
	mkdir ${save_path}/combat_objs
	fi

#GET CSV FILENAME
csv_fname=$(basename $csv .csv)

#######################################################################################
#SET COVAR COLS
covar_list="age_days,sexMale"

#LIST POSSIBLE CONFIGS
config_list="cf.gamlss cf cf.lm cf.gam"
#######################################################################################
for config in $config_list
do
	echo "Prepping $config"
	#write bash script
	bash_script=$bash_dir/${csv_fname}_${config}_combat.sh
	touch $bash_script

	#SIMPLE COMBAT, NO COVARS
	if [ $config = "cf" ]
	then
		echo "singularity run --cleanenv $img Rscript --save $cf_script $csv $batch $save_path $config" > $bash_script

	#COMBAT LM
	elif [ $config = "cf.lm" ]
	then
		echo "singularity run --cleanenv $img Rscript --save $cf_script $csv $batch $save_path $config $covar_list 'lm, formula = y ~ age_days + sexMale'" > $bash_script

	#COMBAT GAM
	elif [ $config = "cf.gam" ]
	then
		echo "singularity run --cleanenv $img Rscript --save $cf_script $csv $batch $save_path $config $covar_list 'gam, formula = y ~ s(age_days) + sexMale'" > $bash_script

	#COMBAT GAMLSS
	elif [ $config = "cf.gamlss" ]
	then
		echo "singularity run --cleanenv $img Rscript --save $cf_script $csv $batch $save_path $config $covar_list 'gamlss, formula = y ~ pb(age_days) + sexMale, sigma.formula = ~ age_days + sexMale'" > $bash_script
	fi

#qsub bash script
	qsub -l h_vmem=64G,s_vmem=64G -N ${config}.${batch}.${csv_fname} -o $bash_dir/${config}.${csv_fname}_${batch}_out.txt -e $bash_dir/${config}.${csv_fname}_${batch}_err.txt $bash_script
done
