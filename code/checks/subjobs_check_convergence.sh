#!/bin/bash
#
#SBATCH --job-name=check_converge
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-482

CONFIGFN=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/config_files/converge_check_config.txt

DF=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIGFN )
PHENO=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $CONFIGFN )
echo "DF: $DF"
echo "PHENO: $PHENO"

# Singularity image path 
SINGULARITY_IMAGE="/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/containers/r_gamlss_0.2.4.sif"

#------------------
BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/

# Run R script inside Singularity container
singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $BASE/code/check_convergence.R $DF $PHENO

# Done!
echo "Job finished running!"
