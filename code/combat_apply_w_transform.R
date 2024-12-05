#script for applying different versions of combat to different feature types
#added optional pipe to qsub gamlss object fitting

#expects 4 - 6 arguments:
## 1. dataframe of data to be combatted
## 2. RDS listing phenotypes to look for in csv
## 3. name of column containing batch identifier, or path to a csv containing batch identifier
## 4. path to save output csv
## 5. list of columns to be included as covariates (OPTIONAL)
## 6. Additional comfam() arguments, including model, formula, ref.batch, ... (OPTIONAL)

set.seed(12345)

#LOAD PACKAGES
library(data.table) 
library(readr)
library(ggplot2) 
library(tidyverse) 
library(mgcv) 
library(gamlss)
library(ComBatFamily)

##########################################################################

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
raw.df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "")
feature_list <- readRDS(args[2])
batch.arg <- args[3]
#if batch arg is csv, merge csv into raw.df and designate last col as batch ID
if (endsWith(batch.arg, '.csv')){
  batch.df <- fread(batch.arg, stringsAsFactors = TRUE, na.strings = "")
  batch <- as.factor(batch.df[,ncol(batch.df)])
  raw.df <- base::merge(raw.df, batch.df)
} else {
  #if batch arg is a column name, select that column name from raw.df
  batch.col <- as.character(batch.arg)
  batch <- as.factor(raw.df[[batch.col]])
}

save_path <- as.character(args[4]) #path to save outputs

#DEF COVARS
covar.list <- as.character(unlist(strsplit(args[5], ",")))
covar.df <- raw.df %>%
  dplyr::select(all_of(covar.list))

stopifnot(length(batch) == nrow(covar.df))
#covar df only works with numeric variables (otherwise need to use matrix to dummy-code). sticking with just df for now, can update later
stopifnot(all(sapply(covar.df, is.numeric)))

cf.args <- args[6]

#extract csv name from input data
csv_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(as.character(args[1])))
csv_basename <- gsub("_", "-", csv_basename)

##########################################################################

#COMBAT

pheno.df <- raw.df %>%
  dplyr::select(any_of(feature_list))
  #replace any 0s w/ 1 pre-log-transform
  pheno.df <- replace(pheno.df, pheno.df==0, 1)
  
  #log-transform ALL pheno vals
  pheno.df <- pheno.df %>%
    mutate(across(c(l), \(x) log(x, base=10)))
  
  #make sure batch, covars, and pheno dfs are all the same length
  stopifnot(nrow(pheno.df) == length(batch))
  
  #turn off empirical bayes if combatting global voluems
  if (grepl("global_vols", args[2]) == TRUE ) {
    eb_arg <- FALSE
  } else {
    eb_arg <- TRUE
  }
  
  #def combatls fun - run CG if RS fails
  cf_gamlss_try <- function(x, eb_arg){
    result <- tryCatch({
      eval(parse(text = paste0("comfam(pheno.df, batch, covar.df, eb = ", eb_arg, ", ", x, ")")))
      } , warning = function(w) {
      message("warning")
      eval(parse(text = paste0("comfam(pheno.df, batch, covar.df, eb = ", eb_arg, ", ", x, ")")))
      } , error = function(e) {
      message("error, trying method=CG()")
      eval(parse(text = paste0("comfam(pheno.df, batch, covar.df, eb = ", eb_arg, ", ", x, ", method = CG())")))
      } , finally = {
      message("done")
      } )
    }
  
  cf.obj <- cf_gamlss_try(cf.args, eb_arg)

  #un-log-transform vals
  cf.obj$dat.combat <- cf.obj$dat.combat %>%
     mutate(across(c(l), \(x) un_log(x)))
  
  #save cf.obj
  saveRDS(cf.obj, file=paste0(save_path, "/combat_objs/", csv_basename, "_cf_obj.rds"))
  
  #row number
  cf.obj.df <- cf.obj$dat.combat %>%
    mutate(id = row_number()) %>%
    dplyr::select(!ends_wit("_X")) #drop non-target cols (used as priors)

#check for negative (impossible) features
total_negative_values <- sum(cf.obj.df < 0)
if (total_negative_values > 0) {
  print(paste("WARNING!", total_negative_values, "negative values found across the following features:"))
  #get names of features with neg. values
  columns_with_negatives <- names(cf.obj.df)[colSums(cf.obj.df < 0) > 0]
  print(columns_with_negatives)
} else {
  print("ComBat successful")
}

#merge back into the rest of the raw dataset (demographics, etc.)
nonpheno.df <- raw.df %>%
  dplyr::select(!any_of(feature_list)) %>%
  mutate(id = row_number())

final.df <- base::merge(cf.obj.df, nonpheno.df, by = "id")

##########################################################################
#WRITE OUT

#append config name
datafile <- paste0(save_path, "/", csv_basename, "_log-cf.gam_batch.", batch.col, "_data.csv")
fwrite(final.df, file=datafile)

print("DONE")
