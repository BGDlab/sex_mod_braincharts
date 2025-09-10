#!/bin/bash
#
#SBATCH --job-name=combatLS
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --array=1
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/combatls_%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/jobfiles/combatls_%A_%a.err

CONFIGFN=$1

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
PHENO=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' "$CONFIGFN")
DF=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
PLIST=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
BATCH=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
COVARS=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )
MMODEL=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $CONFIGFN )
SMODEL=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $CONFIGFN )
SAVE_PATH=$(awk -F'\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $CONFIGFN )

echo "DF: $DF"
echo "PHENO LIST: $PLIST"
echo "BATCH: $BATCH"
echo "COVARS: $COVARS"
echo "MU MODEL: $MMODEL"
echo "SIGMA MODEL: $SMODEL"
echo "SAVE_PATH: $SAVE_PATH"


#------------------

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/

module load R/4.4.0

Rscript $BASE/code/combat_apply_w_transform.R $DF $PLIST $BATCH $COVARS "$MMODEL" "$SMODEL" $SAVE_PATH

# Done!
echo "Job finished running!"
