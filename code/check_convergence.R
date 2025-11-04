#!/usr/bin/env Rscript

#trying to check if models converged successfully
#written with chat gpt

# Load packages
library(dplyr)
library(gamlss)

# Define the directory containing your .RData or .rds files
base_dir_pattern <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/cv_sample_?_train/*/model_objs/"

# Find all directories matching the pattern
dirs <- Sys.glob(base_dir_pattern)

# Find all RDS files ending in "BestMod" within those directories
mod_files <- unlist(lapply(dirs, function(d) {
  list.files(d, pattern = "BestMod\\.rds$", full.names = TRUE)
}))

# Initialize results list
results <- vector("list", length(mod_files))

# Loop through files
for (i in seq_along(mod_files)) {
  file <- mod_files[i]
  
  # Use tryCatch to handle any issues safely
  res <- tryCatch({
    # Load efficiently depending on file type
    obj <- readRDS(file)
    model
    # Extract converged flag (assumes object$converged exists)
    converged_flag <- isTRUE(obj$converged)
    
    cv_sample_dir <- sub(".*(cv_sample[^/]+).*", "\\1", file)
    
    # Return one-row data.frame
    data.frame(
      cv_sample = cv_sample_dir,
      file = basename(file),
      converged = converged_flag,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    # In case of failure
    data.frame(
      cv_sample = cv_sample_dir,
      file = basename(file),
      converged = NA,
      stringsAsFactors = FALSE
    )
  })
  
  results[[i]] <- res
  
  # Clean up memory
  rm(obj)
  gc(verbose = FALSE)
}

# Bind all rows together
df_results <- bind_rows(results)

# Write to CSV
write.csv(df_results, file = file.path("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts", "training_convergence_summary.csv"), row.names = FALSE)
