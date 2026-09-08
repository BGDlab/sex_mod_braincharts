#!/bin/bash
#
#SBATCH --job-name=clean_mods
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/clean_mods/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/clean_mods/R-%A_%a.err


BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts
SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.2.16.sif"
script=$BASE/sharing_code/clean_models_to_share.R
CONFIGFN=$BASE/code/config_files/cv_sample_logAgeTRUE_clean_config.txt

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
MODEL=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
CAT=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
TOTAL=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )
SPLIT=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $CONFIGFN )
SAVEMOD=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $CONFIGFN )
SAVEFIG=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $CONFIGFN )

echo "DF: $DF"
echo "FULL MODEL: $MODEL"
echo "PHENO CATEGORY: $CAT"
echo "TOTAL-SIZE CORRECTED: $TOTAL"
echo "CV SPLIT: $SPLIT"
echo "MODEL SAVE PATH: $SAVEMOD"
echo "FIGURE SAVE PATH: $SAVEFIG"

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $script $DF $MODEL $CAT $TOTAL $SPLIT $SAVEMOD $SAVEFIG

# Done!
echo "Job finished running!"
