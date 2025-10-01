#quick test for adding 1000 to phenos

set.seed(9999)

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

source("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/gamlss_fit_funs.R")

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
pheno <- as.character(args[2])
fs <- as.character(args[3])
total <- as.character(args[4])
save_path <- as.character(args[5])
log_pheno <- as.logical(args[6])
log_age <- as.logical(args[7])
sm <- as.character(args[8])

pheno_sym <- sym(pheno)
  
    
test_mod <- gamlss_lambda(pheno,
                          lambda=NULL,
                          fs_ver = fs,
                          fs_moment="both",
                          fam="BCCGo",
                          weight= FALSE,
                          nu_form="logAge_days",
                          start.from=NULL)
  
  if(is.null(test_mod)){
   stop("Error: model did not converge")
  }
  
  saveRDS(test_mod, paste0("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/test_BCCGo/", pheno, "_test_m2_bccgo_mod.rds"))
  
  #re-write call info to be safe
  test_mod$call$data <- "df"
  test_mod$call$family <- "BCCGo"
  
  fan_plot <- make_centile_fan(test_mod, df, x_var="logAge_days", color_var="sexMale")
  ggsave(file=paste0("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/test_BCCGo/", pheno, "_test_m2_bccgo_centiles.png"), fan_plot)