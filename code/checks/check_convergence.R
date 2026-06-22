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
pheno_val <- as.character(args[2])

# Find all RDS files within those directories
mod_files <- list.files(dir, pattern = paste0("^", pheno_val, ".*\\.rds$"), recursive = TRUE, full.names = TRUE)
n <- length(mod_files)
print(n)

# Initialize results 
df_results <- data.frame(
  cv_sample = character(n),
  file = character(n),
  converged = logical(n),
  y = character(n)
)

# Loop through files
for (i in seq_len(n)) {
  
  obj <- tryCatch({readRDS(file)
  }, error = function(e) NA)
  
  file <- mod_files[i]
  cv_sample_dir <- sub(".*(cv_sample[^/]+).*", "\\1", file)
  df_results$cv_sample[i] <- cv_sample_dir
  df_results$file[i] <- file
  
  df_results$converged[i] <- tryCatch({isTRUE(obj$converged)}, 
                                      error = function(e) NA)
  df_results$y[i] <- tryCatch({
    as.character(obj$mu.terms[[2]])},
    error = function(e) NA_character_)
}

# Write to CSV
df_results$pheno <- pheno_val
csv_name <- paste0(pheno_val, "_convergence_summary.csv")
write.csv(df_results, file = file.path(dir, csv_name), row.names = FALSE)
