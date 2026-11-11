#Select runner-up (2nd-lowest BIC) training model for each pheno, to test w/ fit_cv_mods_test.R as a parallel check on BestMods
#writes manifest read by config_cv_mods_test_runnerup.sh
#run from outside code dir

#LOAD PACKAGES
library(data.table)
library(dplyr)

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
total <- as.character(args[1]) #TRUE/FALSE
log_age <- as.character(args[2]) #TRUE/FALSE
stopifnot(total %in% c("TRUE", "FALSE"), log_age %in% c("TRUE", "FALSE"))

#PATHS
config_path <- "./code/config_files"
pheno_lists <- list.files("./pheno_lists", pattern="\\.txt$", full.names=TRUE)
out_file <- paste0(config_path, "/cv_runnerup_manifest_total", total, "_logAge", log_age, ".csv")

#rebuild model name from summary row - must match m_name in fit_cv_mods.R / fit_cv_total_mods.R
make_m_name <- function(row_df){
  if ("total_moment" %in% names(row_df)){
    paste(paste0("fs", row_df$fs_moment), row_df$nu, paste0("total", row_df$total_moment), sep="_")
  } else {
    paste(row_df$fs_moment, row_df$nu, "BCTo", sep="_")
  }
}

manifest <- data.frame()

#LOOP THROUGH TRAINING SPLITS
for (train_split in c("A", "B")){
  test_split <- ifelse(train_split == "A", "B", "A")
  train_root <- paste0("./cv_sample_", train_split, "_train")

  #LOOP THROUGH PHENO CATEGORIES
  for (pheno_list in pheno_lists){
    pheno_cat <- sub("\\.txt$", "", basename(pheno_list))
    phenos <- readLines(pheno_list)
    phenos <- phenos[phenos != ""]

    #find training dir - allow for old _logPheno*_ naming; anchored so weighted_/age2plus_ dirs are excluded
    dir_pattern <- paste0("^", pheno_cat, "_total", total, "_(.*_)?logAge", log_age, "_pbmods$")
    train_dirs <- list.files(train_root, pattern=dir_pattern, full.names=TRUE)
    if (length(train_dirs) == 0){
      warning(paste("no training dir found in", train_root, "for", pheno_cat))
      next
    } else if (length(train_dirs) > 1){
      stop(paste("multiple training dirs found:", paste(train_dirs, collapse=", ")))
    }
    train_dir <- normalizePath(train_dirs)
    print(paste("searching", train_dir))

    #LOOP THROUGH PHENOS
    for (pheno in phenos){
      sum_file <- file.path(train_dir, "model_sums", paste0(pheno, "_summary.csv"))

      n_mods <- 0L
      best_name <- NA_character_
      ru_name <- NA_character_
      best_bic <- NA_real_
      ru_bic <- NA_real_
      ru_mod <- NA_character_

      if (!file.exists(sum_file)){
        status <- "no_summary"
      } else {
        sum_df <- fread(sum_file) %>%
          arrange(BIC) #same ordering used to pick BestMod in training
        n_mods <- nrow(sum_df)
        best_name <- make_m_name(sum_df[1, ])
        best_bic <- sum_df$BIC[1]

        if (n_mods < 2){
          status <- "lt2_models"
        } else {
          ru_name <- make_m_name(sum_df[2, ])
          ru_bic <- sum_df$BIC[2]
          ru_mod <- file.path(train_dir, "model_objs", paste0(pheno, "_", ru_name, "_mod.rds"))
          best_mod <- file.path(train_dir, "model_objs", paste0(pheno, "_", best_name, "_BestMod.rds"))

          #check summary agrees w/ what's on disk
          if (!file.exists(best_mod)){
            status <- "best_mismatch"
          } else if (!file.exists(ru_mod)){
            status <- "runnerup_missing"
          } else {
            status <- "ok"
          }
        }
      }

      #first 6 cols are read positionally by config_cv_mods_test_runnerup.sh - don't reorder
      manifest <- rbind(manifest, data.frame(
        "train_split" = train_split,
        "test_split" = test_split,
        "pheno_cat" = pheno_cat,
        "pheno" = pheno,
        "status" = status,
        "runnerup_mod" = ru_mod,
        "best_m_name" = best_name,
        "runnerup_m_name" = ru_name,
        "best_BIC" = best_bic,
        "runnerup_BIC" = ru_bic,
        "delta_BIC" = ru_bic - best_bic,
        "n_models" = n_mods,
        "train_dir" = train_dir
      ))
    }
  }
}

stopifnot(nrow(manifest) > 0)
print(table(manifest$train_split, manifest$status))

#SAVE
dir.create(config_path, showWarnings = FALSE)
print(paste("saving to", out_file))
fwrite(manifest, file=out_file)
