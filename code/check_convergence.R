#!/usr/bin/env Rscript

#trying to check if models converged successfully
#written with chat gpt

# Load packages
library(dplyr)
library(gamlss)

# Define the directory containing your .RData or .rds files
args <- commandArgs(trailingOnly = TRUE)
print(args)
dir <- as.character(args[1])
pheno <- as.character(ars[2])

# Find all RDS files within those directories
mod_files <- list.files(dir, pattern = paste0("*/model_objs/", pheno,"*.rds$"), full.names = TRUE)
n <- length(mod_files)

# Initialize results 
df_results <- data.frame(
  cv_sample = character(n),
  pheno = character(n),
  file = character(n),
  converged = logical(n),
  correct_pheno = logical(n),
  stringsAsFactors = FALSE
)

# Loop through files
for (i in seq_len(n)) {
  file <- mod_files[i]
  cv_sample_dir <- sub(".*(cv_sample[^/]+).*", "\\1", file)
  df_results$cv_sample[i] <- cv_sample_dir
  df_results$file[i] <- basename(file)
  
  df_results$converged[i] <- tryCatch({
    obj <- readRDS(file)
    isTRUE(obj$converged)
  }, error = function(e) NA)
  
  df_results$correct_pheno[i] <- tryCatch({
    obj <- readRDS(file)
    pheno==gamlssModel$mu.terms[[2]]
  }, error = function(e) NA)
}

# Write to CSV
csv_name <- paste0(pheno, "_convergence_summary.csv")
write.csv(df_results, file = file.path(dir, csv_name), row.names = FALSE)
