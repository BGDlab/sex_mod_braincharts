#quick test for adding 1000 to phenos

set.seed(9999)

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

source("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/gamlss_fit_funs.R")

df_full <- fread("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/data/v3_CN_cleaned.csv", stringsAsFactors=TRUE)
pheno_list <- readRDS("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/pheno_lists/global_vols.rds")


for (pheno in pheno_list){
  pheno <- as.character(pheno)
  print(pheno)
  
  pheno_sym <- sym(pheno)
  
  df <- df_full %>%
    dplyr::select(all_of(c(pheno, "fs_version_GM", "logAge_days", "sexMale_x_logAge", "sexMale", "study_site"))) %>%
    na.omit() %>%
    trunc_coverage("logAge_days")
    
  test_mod <- gamlss_3lambda(pheno,
                               lambda=NULL, 
                               fs_ver = "fs_version_GM",
                               fs_moment="both",
                               fam="BCCGo", 
                               weight= FALSE,
                               nu_form="logAge_days",
                               start.from=NULL)
  
  if(is.null(test_mod)){
    print("next phenotype")
    next
  }
  
  saveRDS(test_mod, paste0("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/", pheno, "_test_m2_mod.rds"))
  
  #re-write call info to be safe
  test_mod$call$data <- "df"
  test_mod$call$family <- "BCCG"
  
  fan_plot <- make_centile_fan(test_mod, df, x_var="logAge_days", color_var="sexMale")
  ggsave(file=paste0("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/", pheno, "_test_m2_bccgo_centiles.png"), fan_plot)
  
}