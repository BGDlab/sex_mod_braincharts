#!/bin/bash
#
#SBATCH --job-name=cv_mods
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G #almost all jobs run with 8G, upping for resubmission
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/code/jobfiles/cv_train/R-%A_%a.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/code/jobfiles/cv_train/R-%A_%a.err

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc
source $BASE/code/check_branch.sh

CONFIGFN=$BASE/code/config_files/cv_mods_totalFALSE_logAgeTRUE_smpb_config.txt

echo "Config file: $CONFIGFN"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#PARSE CONFIG FILE
DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
PHENO=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
FS=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $CONFIGFN )
TOTAL=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $CONFIGFN )
SAVE_PATH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $CONFIGFN )
LOG_AGE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $CONFIGFN )
SMOOTH=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $CONFIGFN )
N_CORES=${SLURM_CPUS_PER_TASK:-1} #cores the R script fits its model grid on

echo "DF: $DF"
echo "PHENO: $PHENO"
echo "FREESURFER COV: $FS"
echo "TOTAL COV: $TOTAL"
echo "SAVE_PATH: $SAVE_PATH"
echo "LOG-SCALE AGE: $LOG_AGE"
echo "SMOOTH: $SMOOTH"
echo "N CORES: $N_CORES"

#------------------

SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.2.18.sif"

script=$BASE/code/fit_cv_mods.R

echo "SCRIPT: $script"

#fit_cv_mods.R fits its model grid on $SLURM_CPUS_PER_TASK cores (memory is per cpu, so it scales too)

#one BLAS/OpenMP thread per forked worker so the workers don't oversubscribe the allocated cores
#(SINGULARITYENV_ vars get through --cleanenv)
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $script $DF $PHENO $FS $TOTAL $SAVE_PATH $LOG_AGE $SMOOTH $N_CORES

chmod 777 -R $SAVE_PATH

# Done!
echo "Job finished running!"
