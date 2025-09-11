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
f_rh <- as.character(args[2])
name <- as.character(args[3])

# clean string
# Step 1: Replace \" with '
f_rh <- gsub('\\"', "'", f_rh)

# Step 2: Unescape backslashes
f_rh <- gsub('\\\\', '', f_rh)

# Step 3: Remove leading/trailing quotes
f_rh_clean <- gsub('^["\']|["\']$', '', f_rh)

# Create formula
f <- formula(paste(pheno, "~", f_rh_clean))

print(paste("pheno:", pheno))
print("formula:")
print(f)

df_path <- paste0(base, "data/pheno_dfs_totalFALSE/", pheno, "_totalFALSE_logPhenoTRUE_logAgeTRUE.csv")
df <- fread(df_path, stringsAsFactors = TRUE, na.strings = "") #path to csv

#add factor sex variable for by models
df <- df %>%
  mutate(sex = as.factor(ifelse(sexMale == 1, "M", "F")))
#unscale pheno
df[[pheno]] <- unscale(df[[pheno]])

mod <- gamlss2(formula=f, data=df, family ="BCCG")

print("fit, saving mod")

saveRDS(mod, file = paste0(base,"/gamlss2_test_mods/", pheno, "_gamlss2_", name, ".rds"))
