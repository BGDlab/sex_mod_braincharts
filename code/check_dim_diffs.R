# Compare nrow(training csv) vs length(mod$mu.lp) for a list of phenos.
# Usage:
#   Rscript check_dim_diffs.R <phenos.txt> <split> <total> <log_age> [base_path]
# Example:
#   Rscript check_dim_diffs.R failed_phenos.txt B TRUE TRUE

suppressMessages({
  library(data.table); library(dplyr); library(tibble)
})

args     <- commandArgs(trailingOnly = TRUE)
phenos_f <- args[1]
split    <- args[2]
total    <- args[3]                       # "TRUE" or "FALSE"
log_age  <- args[4]                       # "TRUE" or "FALSE"
BASE     <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts"

phenos <- trimws(readLines(phenos_f))
phenos <- phenos[nzchar(phenos)]
cat("Checking", length(phenos), "phenos for split =", split,
    "(total=", total, ", log_age=", log_age, ")\n\n")

dfs_dir   <- file.path(BASE, "data", paste0("cv_sample_", split, "_dfs"))
mods_root <- file.path(BASE, paste0("cv_sample_", split, "_test"))
stopifnot(dir.exists(dfs_dir), dir.exists(mods_root))

# Pre-walk the model tree once (faster than per-pheno globs across 74 phenos)
all_rds <- list.files(mods_root, pattern = "_full_mod\\.rds$",
                      recursive = TRUE, full.names = TRUE)
all_rds <- all_rds[
  grepl(paste0("_total", total, ".*logAge", log_age, "_pbmods/model_objs/"), all_rds) &
    !grepl("weighted", all_rds)
]
cat("Indexed", length(all_rds), "candidate model RDS files\n\n")

results <- lapply(phenos, function(p) {
  csv_pat   <- paste0("^", p, "_total", total, "_logAge", log_age, "\\.csv$")
  csv_match <- list.files(dfs_dir, pattern = csv_pat, full.names = TRUE)
  rds_match <- all_rds[grepl(paste0("/", p, "_[^/]*_full_mod\\.rds$"), all_rds)]
  
  if (length(csv_match) != 1 || length(rds_match) != 1) {
    return(tibble(
      pheno = p, csv_n = NA_integer_, mod_n = NA_integer_, diff = NA_integer_,
      note  = sprintf("csv:%d, rds:%d", length(csv_match), length(rds_match))
    ))
  }
  
  csv_n <- tryCatch(nrow(fread(csv_match, select = 1L)),
                    error = function(e) NA_integer_)
  mod   <- tryCatch(readRDS(rds_match), error = function(e) NULL)
  if (is.null(mod)) {
    return(tibble(pheno = p, csv_n = csv_n, mod_n = NA_integer_,
                  diff = NA_integer_, note = "readRDS failed"))
  }
  tibble(pheno = p, csv_n = csv_n, mod_n = length(mod$mu.lp),
         diff = csv_n - length(mod$mu.lp), note = "")
})

out <- bind_rows(results)

cat("--- Frequency table of (csv_n - mod_n) ---\n")
print(sort(table(out$diff, useNA = "always"), decreasing = TRUE))

cat("\n--- Summary ---\n")
print(summary(out$diff))
cat("Max |diff|:", suppressWarnings(max(abs(out$diff), na.rm = TRUE)), "\n")

cat("\n--- Phenos with diff == 0 (likely a different failure mode) ---\n")
print(out %>% filter(diff == 0) %>% select(pheno))

cat("\n--- Phenos that couldn't be matched ---\n")
print(out %>% filter(is.na(diff)) %>% select(pheno, note))

fwrite(out, "dim_diffs.csv")
cat("\nWrote dim_diffs.csv (", nrow(out), "rows)\n")