#!/bin/bash

#run grab_pngs.R as a job submission
module load R/4.4.0
Rscript /mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/grab_pngs.R $1
