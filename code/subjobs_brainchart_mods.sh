#!/bin/bash
#
#SBATCH --job-name=bc_mods
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G #TEST THIS
#SBATCH --array=1-4#239%30
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/bc_mods/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/bc_mods/R-%A_%a.err

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts
CONFIGFN=$1

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
PHENO=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
MOD_A=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
MOD_B=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $CONFIGFN )

echo "PHENO: $PHENO"
echo "DF: $DF"
echo "MODEL A: $MOD_A"
echo "MODEL B: $MOD_B"
echo "SAVE_PATH: $SAVE_PATH"

#------------------

SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.1.2.sif"

script=$BASE/code/fit_brainchart_mods.R

echo "SCRIPT: $script"

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $script $PHENO $DF $MOD_A $MOD_B $SAVE_PATH

chmod 777 -R $SAVE_PATH

# Done!
echo "Job finished running!"
