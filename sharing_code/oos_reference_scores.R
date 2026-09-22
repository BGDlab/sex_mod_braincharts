# Script for obtaining reference (centile/z scores) from new data 
# benchmarked on brain charts from Gardner et al.

##################################
### LOAD PACKAGES
##################################

library(dplyr)
library(data.table)
library(gamlss)
# install.packages("devtools")
devtools::install_github("BGDlab/gamlssTools@dev", build_vignettes = FALSE) #currently dev version is required
library(gamlssTools)
remotes::install_github("andy1764/gamlss2charts@dev", build_vignettes = FALSE) #currently dev version is required
library(gamlss2charts)
devtools::source_url("https://raw.githubusercontent.com/BGDlab/sex_mod_braincharts/sharing/sharing_code/oos_reference_scores_helper_funs.R") #source helper funs

##################################
### DEFINE ARGUMENTS
##################################

# phenotypes to score - defaults to all
pheno_list <- readLines("https://raw.githubusercontent.com/BGDlab/sex_mod_braincharts/refs/heads/sharing/sharing_code/all_phenos.txt", warn = FALSE)
#subset if desired
# test: 20 random phenotypes
# pheno_list <- sample(pheno_list, 100)

#logical indicating whether to score controlling for total brain size
total <- TRUE

#data to score - set as simulated data template
df <- fread("sharing_code/testdata_fs7/test_output.csv")

batch <- "study_site" # variable containing new levels - batch effects to be estimated and removed

#reference data OR reference condition
ref_data <- ~ dx == "CN"

min_ref <- 75

#filename to save outputs under
out_file <- "brainchart_ref_scores.csv"

#where to get the reference models from.
#  NULL  -> stream each model from GitHub as it is needed (nothing to download,
#           needs a working connection for the whole run)
#  path  -> read models from a local copy of models_to_share/ (download once,
#           then runs offline and is faster); model_ref is ignored
model_dir <- "models_to_share"

# BGDlab/sex_mod_braincharts branch, tag, or commit SHA corresponding to models for reproducibiltiy
# only used when model_dir is NULL
model_ref  <- "main"

##################################
### VALIDATE INPUT DATA
##################################

#check variables used in all phenotype models
required <- c("study_site", "sexMale", "logAge_days")
missing_required <- setdiff(required, names(df))

fs_version_present <- any(c("fs_version_SA", "fs_version_CT", "fs_version_GM") %in% names(df))

msgs <- c(
  if (length(missing_required) > 0)
    paste("Missing required column(s):", paste(missing_required, collapse = ", ")),
  if (!fs_version_present)
    "None of fs_version_SA, fs_version_CT, or fs_version_GM found in df"
)

if (length(msgs) > 0) stop(paste(msgs, collapse = "\n"))

stopifnot(
  "sexMale must be numeric" = is.numeric(df$sexMale),
  "sexMale must contain only 0/1" = all(df$sexMale %in% c(0, 1))
)

df <- df %>%
  mutate(sexMale_x_logAge = sexMale * logAge_days, #calculate sex x age interaction
         .row_id = seq_len(n())) #stable row key

#with a local model directory, fail now if it is missing or incomplete
check_model_dir(pheno_list, total)

##################################
### CALCULATE REFERENCE SCORES
##################################
print("calculating centiles...")

df_cent <- Filter(Negate(is.null),
                  lapply(pheno_list, score_pheno, df = df, 
                         total = total, 
                         batch = batch, 
                         ref_data = ref_data, min_ref = min_ref))

#rejoin to the full input data by row key
df_full_cent <- Reduce(
  function(a, b) dplyr::left_join(a, b, by = ".row_id"),
  df_cent,
  init = as.data.frame(df)
) %>%
  dplyr::select(-.row_id)

print(paste0("scored ", length(df_cent), "/", length(pheno_list), " phenotypes on ",
             nrow(df_full_cent), " subjects"))

##################################
### SAVE OUTPUTS
##################################
fwrite(df_full_cent, out_file)
