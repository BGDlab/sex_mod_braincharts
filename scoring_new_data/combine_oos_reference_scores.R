# Stitch per-phenotype outputs from submit_oos_reference_scores.sh into one file
#   Rscript combine_oos_reference_scores.R <work_dir> <out_file>

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript combine_oos_reference_scores.R <work_dir> <out_file>")
work_dir <- args[1]
out_file <- args[2]

phenos <- readLines(file.path(work_dir, "phenos.txt"), warn = FALSE)
files <- file.path(work_dir, "scores", paste0(phenos, ".csv"))

have <- file.exists(files)
if (!any(have)) stop("no per-phenotype outputs found in ", file.path(work_dir, "scores"))
if (any(!have))
  warning(sum(!have), " phenotype job(s) wrote no output (check logs in ",
          file.path(work_dir, "logs"), "): ", paste(phenos[!have], collapse = ", "),
          call. = FALSE)
phenos <- phenos[have]
files <- files[have]

score_cols <- function(p) paste0(p, c("_z", "_centile"))

#start from the copy of --df every task scored, whose rows were numbered (.row_id)
#by submit_oos_reference_scores.sh
df_copy <- file.path(work_dir, "df.csv")
if (!file.exists(df_copy)) stop("copy of --df not found: ", df_copy, call. = FALSE)
out <- fread(df_copy, na.strings = c("NA", "", '""'))
if (!".row_id" %in% names(out)) stop("no .row_id column in ", df_copy, call. = FALSE)
n <- nrow(out)

#each per-pheno file has .row_id + that pheno's score columns for the rows it scored;
#match the scores to rows by .row_id, leaving NA for rows that weren't scored
for (i in seq_along(files)) {
  s <- fread(files[i])
  cols <- intersect(score_cols(phenos[i]), names(s))
  if (length(cols) == 0 || nrow(s) == 0) {
    message(phenos[i], ": no score columns (scoring failed; see its log)")
    next
  }
  if (!".row_id" %in% names(s))
    stop(phenos[i], ": no .row_id column in ", files[i], call. = FALSE)
  rows <- match(s$.row_id, out$.row_id)
  if (anyNA(rows))
    stop(phenos[i], ": ", sum(is.na(rows)), " .row_id value(s) not in ", df_copy, call. = FALSE)
  if (anyDuplicated(rows))
    stop(phenos[i], ": duplicated .row_id values", call. = FALSE)

  set(out, j = cols, value = NA_real_)
  set(out, i = rows, j = cols, value = s[, cols, with = FALSE])
}

out[, .row_id := NULL]
fwrite(out, out_file)
n_scored <- sum(paste0(phenos, "_centile") %in% names(out))
cat("combined ", n_scored, "/", length(have), " phenotypes on ", n,
    " subjects -> ", out_file, "\n", sep = "")
