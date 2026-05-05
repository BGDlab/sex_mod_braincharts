#Check that every total-corrected BestMod actually used the total-size covariate
#in at least one moment (i.e., total_moment != "none").
#
#Usage: Rscript code/check_total_in_bestmods.R [base_dir]
#       base_dir defaults to the current working directory.

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else "."

pattern <- file.path(
  base_dir,
  "cv_sample_*_train",
  "*_totalTRUE_logAge*_*mods",
  "model_objs",
  "*_BestMod.rds"
)
files <- Sys.glob(pattern)

if (length(files) == 0) {
  stop("No total-corrected BestMod files found under: ", base_dir)
}

total_moment <- sub(".*_total([a-z]+)_BestMod\\.rds$", "\\1", basename(files))
none <- files[total_moment == "none"]

cat(sprintf("Found %d total-corrected BestMod files\n", length(files)))
cat("total_moment counts:\n"); print(table(total_moment, useNA = "ifany"))

if (length(none) == 0) {
  cat("\nAll BestMods include the total-size covariate in at least one moment.\n")
} else {
  cat(sprintf("\n%d BestMod(s) had total_moment == 'none':\n", length(none)))
  cat(paste0("  ", none, collapse = "\n"), "\n")
}
