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
devtools::source_url("https://githubusercontent.com") #source helper funs

##################################
### DEFINE ARGUMENTS
##################################

# phenotypes to score - any subset of phenos modeled, defaults to all
full_pheno_list <- "https://raw.githubusercontent.com/BGDlab/sex_mod_braincharts/refs/heads/main/sharing_code/all_phenos.txt"
pheno_list <- readLines(full_pheno_list, warn = FALSE) %>% as.list()

#logical indicating whether to score controlling for total brain size
total <- FALSE

#data to score - set as simulated data template
df <- fread()

batch <- "study_site" # variable containing new levels - batch effects to be estimated and removed

#reference data OR reference condition

#filename to save outputs under

# BGDlab/sex_mod_braincharts branch, tag, or commit SHA corresponding to models for reproducibiltiy
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

##################################
### CALCULATE REFERENCE SCORES
##################################
print("calculating centiles...")

df_cent <- Filter(Negate(is.null),
                  lapply(pheno_list, score_pheno, df = df, total = total, batch = batch))

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
