#seeing if models behave better when age is not log-scaled
set.seed(99999)

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

source("./code/gamlss_fit_funs.R")

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
pheno <- as.character(args[2]) #phenotype
fs <- as.character(args[3]) #freesurfer version
total <- as.character(args[4]) # total pheno to control for (or NULL)
log_age <- as.logical(args[5]) #log-scale age TRUE/FALSE
save_path <- as.character(args[6]) #path to save df

#make completely sure NDAR-SCZGlu has been removed
df <- df %>%
  filter(study != "NDAR-SCZGlu" & study != "NDAR_SCZGlu")

#log-scale age if necessary
if (log_age == TRUE){
  age_var <- "logAge_days"
  sex_age_var <- "sexMale_x_logAge"
} else if (log_age == FALSE) {
  age_var <- "age_days"
  sex_age_var <- "sexMale_x_age"
}

#drop extra variables
if (total == 'NULL'){
  df <- df %>%
    dplyr::select(any_of(c(pheno, fs, age_var, sex_age_var, "sexMale", "study_site"))) %>%
    na.omit() %>%
    trunc_coverage(age_var, max_loops=100) #drop points at ends if too sparse
} else {
  df <- df %>%
    dplyr::select(any_of(c(pheno, fs, age_var, sex_age_var, "sexMale", "study_site", total))) %>%
    na.omit() %>%
    trunc_coverage(c(total, age_var), max_loops=100) #drop points at ends if too sparse
}

print(names(df))

#total name
if (total=="NULL"){
  total_tf <- "FALSE"
} else {
  total_tf <- "TRUE"
}


print(unique(df$study))

fwrite(df, file=paste0(save_path, "/", pheno, "_total", total_tf, "_logAge" , log_age, ".csv"))

print("DONE")