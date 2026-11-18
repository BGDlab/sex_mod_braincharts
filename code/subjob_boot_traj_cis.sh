#!/bin/bash
#
#SBATCH --job-name=boot_traj
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/code/jobfiles/combine_boot/R-traj-%A.out
#SBATCH --error=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/code/jobfiles/combine_boot/R-traj-%A.err
#usage: sbatch ./code/subjob_boot_traj_cis.sh <split A boot dir> <split B boot dir> <out .rds>
#e.g. sbatch ./code/subjob_boot_traj_cis.sh ./cv_sample_A_boot/global_vols_totalTRUE_logAgeTRUE_boot \
#       ./cv_sample_B_boot/global_vols_totalTRUE_logAgeTRUE_boot ./data/glob_traj_boot_cis.rds

BOOT_DIR_A=$(realpath $1)
BOOT_DIR_B=$(realpath $2)
OUT_FILE=$(realpath -m $3)

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc
source $BASE/code/check_branch.sh

SINGULARITY_IMAGE="$BASE/containers/r_gamlss_0.2.19.sif"

script=$BASE/code/boot_traj_cis.R

echo "SCRIPT: $script"
echo "BOOT_DIR_A: $BOOT_DIR_A"
echo "BOOT_DIR_B: $BOOT_DIR_B"
echo "OUT_FILE: $OUT_FILE"

#output readable/writable by everyone, whoever submits the job
umask 000

singularity run --cleanenv \
    -B $BASE \
    $SINGULARITY_IMAGE \
    Rscript $script $BOOT_DIR_A $BOOT_DIR_B $OUT_FILE

# Done!
echo "Job finished running!"
