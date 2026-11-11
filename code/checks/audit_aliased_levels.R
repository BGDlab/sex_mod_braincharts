#Audit shared models for aliased (NA) fixed-effect coefficients caused by empty factor levels
#(e.g. fs_version level "" from quoted empty strings in the input csv that survived as a factor
#level after na.omit() removed all its rows). An empty reference level makes the remaining level
#dummies sum to the intercept, so the last level's coefficient is NA. predict.gamlss() tolerates
#this, but out-of-sample scoring of new batches (gamlss2charts data-free offsets) returns NA for
#every row.
#
#Usage: Rscript code/checks/audit_aliased_levels.R <model_dir> [out_csv]

library(gamlss)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript audit_aliased_levels.R <model_dir> [out_csv]")
model_dir <- args[1]
out_csv <- if (length(args) >= 2) args[2] else file.path(model_dir, "aliased_levels_audit.csv")

files <- list.files(model_dir, pattern = "_sharing_model\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("no *_sharing_model.rds files found in ", model_dir)
print(paste("auditing", length(files), "models in", model_dir))

#audit one model: one row per parameter
audit_model <- function(f) {
  m <- tryCatch(readRDS(f), error = function(e) NULL)
  fname <- basename(f)
  meta <- regmatches(fname, regexec("^(.*)_split([^_]+)_total(TRUE|FALSE)_sharing_model\\.rds$", fname))[[1]]
  if (is.null(m)) {
    return(data.frame(file = fname, pheno = meta[2], split = meta[3], total = meta[4],
                      param = NA, na_coefs = NA, empty_levels = NA, ref_levels = NA,
                      read_error = TRUE))
  }

  rows <- lapply(m$parameters, function(p) {
    cf <- coef(m, p)
    #smoother terms (pb/random) have no linear coefficient by design; ignore them
    sm <- colnames(m[[paste0(p, ".s")]])
    na_fixed <- names(cf)[is.na(cf) & !names(cf) %in% sm]

    xlev <- m[[paste0(p, ".xlevels")]]
    empty <- unlist(lapply(names(xlev), function(v)
      if (any(!nzchar(xlev[[v]]))) v))
    ref <- vapply(names(xlev), function(v) paste0(v, "=", xlev[[v]][1]), character(1))

    data.frame(file = fname, pheno = meta[2], split = meta[3], total = meta[4],
               param = p,
               na_coefs = paste(na_fixed, collapse = "; "),
               empty_levels = paste(empty, collapse = "; "), #factors with a "" level
               ref_levels = paste(ref, collapse = "; "),     #first (reference) level of each factor
               read_error = FALSE)
  })
  rm(m); gc(verbose = FALSE)
  do.call(rbind, rows)
}

res <- do.call(rbind, lapply(seq_along(files), function(i) {
  if (i %% 25 == 0) print(paste(i, "/", length(files)))
  audit_model(files[i])
}))

res$aliased <- !is.na(res$na_coefs) & nzchar(res$na_coefs)
res$has_empty_level <- !is.na(res$empty_levels) & nzchar(res$empty_levels)

data.table::fwrite(res, out_csv)

##### SUMMARY #####
by_model <- aggregate(cbind(aliased, has_empty_level, read_error) ~ file + pheno + split + total,
                      data = res, FUN = any)

cat("\n==== SUMMARY ====\n")
cat("models audited:          ", nrow(by_model), "\n")
cat("unreadable:              ", sum(by_model$read_error), "\n")
cat("with NA fixed-effect coef:", sum(by_model$aliased), "\n")
cat("with a \"\" factor level:  ", sum(by_model$has_empty_level), "\n")
cat("NA coef but no \"\" level: ", sum(by_model$aliased & !by_model$has_empty_level),
    "(aliased for some other reason - inspect these)\n")

cat("\naffected models by total/split:\n")
print(with(by_model, table(total = total, split = split, aliased = aliased)))

cat("\nwhich coefficients are NA (count of model x parameter rows):\n")
na_terms <- unlist(strsplit(res$na_coefs[res$aliased], "; "))
print(sort(table(na_terms), decreasing = TRUE))

cat("\nfull results written to", out_csv, "\n")
