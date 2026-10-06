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

#every output is the full input df + that pheno's score columns, in the same row order,
#so keep the input columns once and add each pheno's score columns
first <- fread(files[1], na.strings = c("NA", "", '""'))
out <- first[, setdiff(names(first), score_cols(phenos[1])), with = FALSE]

for (i in seq_along(files)) {
  hdr <- names(fread(files[i], nrows = 0))
  cols <- intersect(score_cols(phenos[i]), hdr)
  if (length(cols) == 0) {
    message(phenos[i], ": no score columns (scoring failed; see its log)")
    next
  }
  s <- fread(files[i], select = cols)
  if (nrow(s) != nrow(out))
    stop(phenos[i], ": ", nrow(s), " rows but expected ", nrow(out), call. = FALSE)
  out[, (cols) := s]
}

fwrite(out, out_file)
n_scored <- sum(paste0(phenos, "_centile") %in% names(out))
cat("combined ", n_scored, "/", length(have), " phenotypes on ", nrow(out),
    " subjects -> ", out_file, "\n", sep = "")
