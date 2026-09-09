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

##################################
### DEFINE ARGUMENTS
##################################

#phenotypes to score - any subset of 
# defaults to all
full_pheno_list <- "https://githubusercontent.com" #need to update
pheno_list <- readLines(full_pheno_list, warn = FALSE) %>% as.list()

#logical indicating whether to score controlling for total brain size
total <- FALSE

#data to score - set as simulated data template
df <- fread()

#reference data OR reference condition

#filename to save outputs under


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

#calculate sex x age interaction
df <- df %>%
  mutate(sexMale_x_logAge = sexMale * logAge_days)

##################################
### CALCULATE REFERENCE SCORES
##################################
print("calculating centiles...")

df_cent <- lapply(pheno_list, function(pheno){
  m <- mod_list[[mn]]
  out_df <- pred_og_centile(
    m,
    og.data = df.og,
    new.data = df_clean,
    get.std.scores = TRUE
  )
  # rename columns dynamically to include pheno and model type
  names(out_df) <- paste0(
    pheno, "_",
    names(out_df), "_",
    mn
  )
  out_df
})


##################################
### SAVE OUTPUTS
##################################