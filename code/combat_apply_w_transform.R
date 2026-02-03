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
library(dplyr) 
library(gamlss)
library(ComBatFamily)
library(splines)
library(gamlssTools)

##########################################################################

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
raw.df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "")
feature_list <- readRDS(args[2])
feature_list <- c(feature_list, paste0(feature_list, "_X"))
batch.col <- as.character(args[3])
#DEF COVARS
covar.list <- as.character(unlist(strsplit(args[4], ",")))

#FILTER DF
filt.df <- raw.df %>%
  dplyr::select(any_of(c(covar.list, feature_list, batch.col, "INDEX.ID"))) %>%
  na.omit() %>%
  group_by(!!sym(batch.col)) %>%
  filter(n() >=5) %>% #remove sites with < 5 ppl
  ungroup()

covar.df <- filt.df %>%
  dplyr::select(all_of(covar.list))

batch <- as.factor(filt.df[[batch.col]])

stopifnot(length(batch) == nrow(covar.df))
#covar df only works with numeric variables (otherwise need to use matrix to dummy-code). sticking with just df for now, can update later
stopifnot(all(sapply(covar.df, is.numeric)))

mu.form <- args[5]
sig.form <- args[6]
# print(mu.form)
# print(sig.form)

save_path <- as.character(args[7]) #path to save outputs

#extract csv name from input data
csv_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(as.character(args[1])))
csv_basename <- gsub("_", "-", csv_basename)

##########################################################################

#COMBAT

df <- filt.df %>%
  dplyr::select(any_of(feature_list))
stopifnot(all(sapply(df, is.numeric)))

  #replace any 0s w/ 1 pre-log-transform
  df <- replace(df, df==0, 1)
  
  #log-transform ALL pheno vals
  df <- df %>%
    mutate(across(any_of(feature_list), \(x) log(x, base=10)))
  
  #make sure batch, covars, and pheno dfs are all the same length
  stopifnot(nrow(df) == length(batch))
  
  #turn off empirical bayes if combatting global voluems
  if (grepl("global_vols", args[2]) == TRUE ) {
    print("gobal volumes, no EB")
    eb_arg <- FALSE
  } else {
    eb_arg <- TRUE
  }
  
  cf.obj <- eval(parse(text = paste0("comfam(df, batch, covar.df, gamlss, formula= y", mu.form, 
                                     ", sigma.formula=", sig.form, 
                                     ", eb = ", eb_arg, 
                                     ", control=gamlss.control(n.cyc=800))")))

  #un-log-transform vals
  cf.obj$dat.combat <- cf.obj$dat.combat %>%
    as.data.frame() %>%
    mutate(across(any_of(feature_list), \(x) un_log(x)))
  
  #save cf.obj
  saveRDS(cf.obj, file=paste0(save_path, "/combat_objs/", csv_basename, "_cf_obj.rds"))
  
  #row number
  cf.obj.df <- cf.obj$dat.combat %>%
    as.data.frame() %>%
    mutate(id = row_number()) %>%
    dplyr::select(!ends_with("_X")) #drop non-target cols (used as priors)

#check for negative (impossible) features
total_negative_values <- sum(cf.obj.df < 0)
total_missing_values <- sum(is.na(cf.obj.df))
if (total_missing_values > 0) {
  print(paste("ERROR!", total_missing_values, "NA values found across the following features:"))
  #get names of features with neg. values
  columns_with_na <- names(cf.obj.df)[colSums(is.na(cf.obj.df)) > 0]
  print(columns_with_na)
} else if (total_negative_values > 0) {
  print(paste("WARNING!", total_negative_values, "negative values found across the following features:"))
  #get names of features with neg. values
  columns_with_negatives <- names(cf.obj.df)[colSums(cf.obj.df < 0) > 0]
  print(columns_with_negatives)
} else {
  print("ComBat successful")
}

#merge back into the rest of the raw dataset (demographics, etc.)
nonpheno.df <- filt.df %>%
  dplyr::select(!any_of(feature_list)) %>%
  mutate(id = row_number())

final.df <- base::merge(cf.obj.df, nonpheno.df, by = "id") %>%
  select(!id) #drop id col post-merge

##########################################################################
#WRITE OUT

#append config name
datafile <- paste0(save_path, "/", csv_basename, "_log-scale_batch.", batch.col, "_data.csv")
fwrite(final.df, file=datafile)

print("DONE")
