#Check already-fit test models against the abs() lambda check now in fit_cv_mods_test*.R
#(old one-sided check only caught refit lambdas smaller than the training lambdas)
#compares each test *_full_mod.rds to the training model it was refit from
#covers main, weighted, age2plus & runner-up test models
#run from outside code dir: Rscript code/checks/check_test_lambdas.R [out_csv]

#LOAD PACKAGES
library(data.table)
library(dplyr)

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
out_file <- ifelse(length(args) > 0, args[1], "./code/checks/test_lambda_check.csv")
tol <- 0.01 #same threshold as fit_cv_mods_test*.R

test_files <- Sys.glob("./cv_sample_?_test*/*/model_objs/*_full_mod.rds")
print(paste(length(test_files), "test models found"))
stopifnot(length(test_files) > 0)

#find training model a test model was refit from - reverses naming in config_cv_mods_test*.sh & fit_cv_mods_test*.R
find_train_mod <- function(test_file){
  test_dir <- basename(dirname(dirname(test_file))) #e.g. global_vols_totalFALSE_logAgeTRUE_pbmods
  test_root <- basename(dirname(dirname(dirname(test_file)))) #e.g. cv_sample_A_test
  test_split <- sub("^cv_sample_([AB])_test.*$", "\\1", test_root)
  train_split <- ifelse(test_split == "A", "B", "A")

  #training dirs may have old _logPheno*_ naming
  dir_pattern <- sub("^(.*_total(TRUE|FALSE))_(logAge(TRUE|FALSE)_pbmods)$", "^\\1_(.*_)?\\3$", test_dir)
  train_dirs <- list.files(paste0("./cv_sample_", train_split, "_train"), pattern=dir_pattern, full.names=TRUE)

  stem <- sub("_full_mod\\.rds$", "", basename(test_file))
  if (grepl("_test_(weighted|age2plus)$", stem)){
    train_name <- sub("_test_(weighted|age2plus)$", "_train_\\1", stem)
  } else if (grepl("_test_runnerup$", test_root)){
    train_name <- sub("_test$", "_mod", stem)
  } else {
    train_name <- sub("_test$", "_BestMod", stem)
  }

  train_files <- file.path(train_dirs, "model_objs", paste0(train_name, ".rds"))
  train_files[file.exists(train_files)]
}

#max abs lambda difference for one moment; Inf if number of smooth terms differs
lambda_diff <- function(base, refit){
  if (length(base) != length(refit)) return(Inf)
  max(abs(base - refit))
}

results <- data.frame()

#LOOP THROUGH TEST MODELS
for (test_file in test_files){
  print(test_file)
  train_file <- find_train_mod(test_file)

  row <- data.frame(
    "test_root" = basename(dirname(dirname(dirname(test_file)))),
    "test_dir" = basename(dirname(dirname(test_file))),
    "test_mod" = basename(test_file),
    "train_mod" = NA_character_,
    "status" = NA_character_,
    "mu_diff" = NA_real_,
    "sigma_diff" = NA_real_,
    "nu_diff" = NA_real_,
    "old_check_pass" = NA
  )

  if (length(train_file) != 1){
    row$status <- ifelse(length(train_file) == 0, "no_train_mod", "multiple_train_mods")
    results <- rbind(results, row)
    next
  }
  row$train_mod <- normalizePath(train_file)

  base_mod <- tryCatch(readRDS(train_file), error = function(e) NULL)
  model <- tryCatch(readRDS(test_file), error = function(e) NULL)
  if (is.null(base_mod) | is.null(model)){
    row$status <- "read_error"
    results <- rbind(results, row)
    next
  }

  row$mu_diff <- lambda_diff(base_mod$mu.lambda, model$mu.lambda)
  row$sigma_diff <- lambda_diff(base_mod$sigma.lambda, model$sigma.lambda)
  signed <- c(base_mod$mu.lambda - model$mu.lambda, base_mod$sigma.lambda - model$sigma.lambda)
  #nu only checked if training model has a nu smooth, as in fit_cv_mods_test*.R
  if (!is.null(base_mod$nu.lambda)){
    row$nu_diff <- lambda_diff(base_mod$nu.lambda, model$nu.lambda)
    signed <- c(signed, base_mod$nu.lambda - model$nu.lambda)
  }

  row$old_check_pass <- all(signed < tol)
  row$status <- ifelse(all(c(row$mu_diff, row$sigma_diff, row$nu_diff) < tol, na.rm=TRUE), "pass", "fail")
  results <- rbind(results, row)

  #get back memory
  rm(base_mod, model)
  gc()
}

#SUMMARY
print(table(results$test_root, results$status))
fails <- results %>% filter(status != "pass")
if (nrow(fails) > 0){
  print("models not passing:")
  print(fails %>% select(test_root, test_dir, test_mod, status, mu_diff, sigma_diff, nu_diff))
} else {
  print("all test models pass")
}

print(paste("saving to", out_file))
fwrite(results, file=out_file)
