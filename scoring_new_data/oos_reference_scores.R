# Script for obtaining reference (centile/z scores) from new data 
# benchmarked on brain charts from Gardner et al.

##################################
### LOAD PACKAGES
##################################

library(dplyr)
library(data.table)
library(gamlss)
if (!requireNamespace("gamlssTools", quietly = TRUE))
  remotes::install_github("BGDlab/gamlssTools@dev", upgrade = "never", build_vignettes = FALSE)
if (!requireNamespace("gamlss2charts", quietly = TRUE))
  remotes::install_github("andy1764/gamlss2charts@dev", upgrade = "never", build_vignettes = FALSE)
library(gamlssTools)
library(gamlss2charts)

##################################
### LOAD HELPER FUNCTIONS
##################################

#source helpers from github; if that fails (e.g. no internet), fall back to
#a local copy next to this script
helper_url <- "https://raw.githubusercontent.com/BGDlab/sex_mod_braincharts/refs/heads/sharing/scoring_new_data/oos_reference_scores_helper_funs.R"

script_dir <- function() {
  #Rscript
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) > 0) return(dirname(normalizePath(gsub("~+~", " ", f[1], fixed = TRUE))))
  #source()'d interactively
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  getwd()
}

tryCatch(
  source(helper_url),
  error = function(e) {
    local_helpers <- file.path(script_dir(), "oos_reference_scores_helper_funs.R")
    if (!file.exists(local_helpers))
      stop("could not source helper functions from GitHub (", conditionMessage(e),
           ") and no local copy found at ", local_helpers, call. = FALSE)
    warning("could not source helper functions from GitHub (", conditionMessage(e),
            "); using local copy at ", local_helpers, call. = FALSE)
    source(local_helpers)
  }
)

##################################
### DEFINE ARGUMENTS
##################################

args <- parse_args(commandArgs(trailingOnly = TRUE))
print(args)

# phenotypes to score - defaults to all
pheno_list <- readLines(args$pheno_list, warn = FALSE)

#data to score - set as simulated data template
df <- fread(args$df)

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
check_model_dir(args$model_dir, pheno_list, args$total)

##################################
### CALCULATE REFERENCE SCORES
##################################
print("calculating centiles...")

df_cent <- Filter(Negate(is.null),
                  lapply(pheno_list, score_pheno, 
                         df = df, 
                         total = args$total, 
                         batch = args$batch, 
                         ref_data = args$ref_data, 
                         min_ref = args$min_ref,
                         model_dir = args$model_dir,
                         model_ref = args$model_ref))

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
fwrite(df_full_cent, args$out_file)
