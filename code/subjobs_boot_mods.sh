#!/bin/bash
#
#SBATCH --job-name=boot_mods
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/code/jobfiles/boot/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/code/jobfiles/boot/R-%A_%a.err

CONFIGFN=$1
LINE=$2

echo "Config file: $CONFIGFN"
echo "Config line: $LINE"
echo "SLURM_ARRAY_TASK_ID (boot id): $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
DF=$(awk -F'\t' -v L=$LINE '$1+0==L {print $2}' $CONFIGFN )
TRAIN_MOD=$(awk -F'\t' -v L=$LINE '$1+0==L {print $3}' $CONFIGFN )
TEST_MOD=$(awk -F'\t' -v L=$LINE '$1+0==L {print $4}' $CONFIGFN )
SAVE_PATH=$(awk -F'\t' -v L=$LINE '$1+0==L {print $5}' $CONFIGFN )
SIMDF=$(awk -F'\t' -v L=$LINE '$1+0==L {print $6}' $CONFIGFN )
SEED=$(awk -F'\t' -v L=$LINE '$1+0==L {print $7}' $CONFIGFN )
STRATA=$(awk -F'\t' -v L=$LINE '$1+0==L {print $8}' $CONFIGFN )

echo "DF: $DF"
echo "TRAINING MODEL: $TRAIN_MOD"
echo "TEST MODEL: $TEST_MOD"
echo "SAVE_PATH: $SAVE_PATH"
echo "SIM DF: $SIMDF"
echo "SEED: $SEED"
echo "STRATA: $STRATA"

#------------------

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/
source $BASE/code/check_branch.sh
#same container as the test fits (fit_cv_mods_test.R)
SINGULARITY_IMAGE="/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/containers/r_gamlss_0.2.18.sif"

#outputs readable/writable by everyone, whoever submits the job
umask 000

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $BASE/code/fit_boot_mods.R $DF $TRAIN_MOD $TEST_MOD $SAVE_PATH $SIMDF $SLURM_ARRAY_TASK_ID $SEED $STRATA

# Done!
echo "Job finished running!"
