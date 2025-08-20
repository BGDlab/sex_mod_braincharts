#quick test for gamlss 2

set.seed(9999)

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlss2)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
pheno <- as.character(args[1])
f_rh <- as.character(args[2]) #formula

f <- formula(paste(pheno, "~", f_rh))

print(paste("pheno:", pheno))
print("formula:")
print(f)

df_path <- paste0(base, "data/pheno_dfs_totalFALSE/", pheno, "_totalFALSE_logPhenoTRUE_logAgeTRUE.csv")
df <- fread(df_path, stringsAsFactors = TRUE, na.strings = "") #path to csv

mod <- gamlss2(formula=f, data=df, family ="BCCG")

print("fit, saving mod")

saveRDS(mod, file = paste0(base, pheno, "_gamlss2_testmod.rds"))
