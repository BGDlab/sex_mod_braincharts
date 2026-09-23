#helper functions for getting out-of-sample reference scores

#argument defaults
DEFAULTS <- list(
  pheno_list = "https://raw.githubusercontent.com/BGDlab/sex_mod_braincharts/refs/heads/sharing/sharing_code/all_phenos.txt", 
  total = NULL,
  df = NULL,
  batch = "study_site",
  ref_data = NULL,
  min_ref = 75,
  out_file = NULL,
  model_dir = NULL, #NULL=stream from github, otherwise local path to models
  model_ref = "sharing" # BGDlab/sex_mod_braincharts branch, tag, or commit SHA corresponding to models for reproducibiltiy
                     # only used when model_dir is NULL
)

# printed by --help
USAGE <- c(
  "Score new data against the Gardner et al. brain charts (z-scores and centiles).",
  "",
  "  Rscript sharing_code/oos_reference_scores.R \\",
  "    --df my_datadict.csv \\",
  "    --total TRUE \\",
  "    --out_file my_ref_scores.csv \\",
  "    [--ref_data \"dx == 'CN'\"] [--min_ref 75] [--batch study_site] \\",
  "    [--pheno_list phenos.txt] [--model_dir models_to_share] [--model_ref main]",
  "",
  "required: --df, --total (TRUE/FALSE), --out_file",
  "--ref_data is either a condition on columns of --df (e.g. \"dx == 'CN'\")",
  "  or a path to a CSV of reference rows",
  "--model_dir reads models from a local copy of models_to_share/;",
  "  leave unset to stream them from GitHub at --model_ref"
)

#parse args
parse_args <- function(argv) {
  #get help
  if (length(argv) == 0 || any(argv %in% c("-h", "--help"))) {
    cat(USAGE, sep = "\n")
    cat("\noptions:\n  --",
        paste(names(DEFAULTS), collapse = "\n  --"), "\n", sep = "")
    # don't take an interactive session down with us
    if (interactive()) stop("no arguments given; see the usage above", call. = FALSE)
    quit(status = if (length(argv) == 0) 1 else 0)
  }
  
  #parse named args
  opts <- DEFAULTS
  seen <- character(0)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
    has_eq <- grepl("=", a, fixed = TRUE)
    flag <- sub("^--", "", if (has_eq) sub("=.*$", "", a) else a)
    key <- gsub("-", "_", flag)   # --min-ref and --min_ref both work
    if (!key %in% names(DEFAULTS)) stop("unknown option: --", flag, call. = FALSE)
    if (key %in% seen) stop("--", flag, " given more than once", call. = FALSE)
    seen <- c(seen, key)
    if (has_eq) {
      val <- sub("^[^=]*=", "", a)
    } else {
      #the next token is the value, unless it is the next flag
      if (i == length(argv) || startsWith(argv[i + 1], "--"))
        stop("missing value for --", flag, call. = FALSE)
      val <- argv[i + 1]; i <- i + 1
    }
    if (!nzchar(val)) stop("empty value for --", flag, call. = FALSE)
    opts[[key]] <- val
    i <- i + 1
  }
  
  #validate parsed args & convert to correct types
  missing <- c("df", "total", "out_file")[vapply(opts[c("df", "total", "out_file")], is.null, logical(1))]
  if (length(missing) > 0)
    stop("missing required option(s): ", paste0("--", missing, collapse = ", "), call. = FALSE)
  
  #total only goes into model filenames, so keep it as the "TRUE"/"FALSE" string
  total <- toupper(opts$total)
  if (!total %in% c("TRUE", "FALSE", "T", "F"))
    stop("--total must be TRUE or FALSE, not '", opts$total, "'", call. = FALSE)
  opts$total <- if (total %in% c("TRUE", "T")) "TRUE" else "FALSE"
  
  #min_ref must be numeric
  min_ref <- suppressWarnings(as.numeric(opts$min_ref))
  if (length(min_ref) != 1 || is.na(min_ref) || min_ref < 0)
    stop("--min_ref must be a non-negative number, not '", opts$min_ref, "'", call. = FALSE)
  opts$min_ref <- min_ref
  
  #files must exist locally or on github
  if (!is_url(opts$df) && !file.exists(opts$df))
    stop("--df file not found: ", opts$df, call. = FALSE)
  if (!is_url(opts$pheno_list) && !file.exists(opts$pheno_list))
    stop("--pheno_list file not found: ", opts$pheno_list, call. = FALSE)
  out_parent <- dirname(opts$out_file)
  if (!dir.exists(out_parent))
    stop("directory for --out_file does not exist: ", out_parent, call. = FALSE)
  
  #ref_data: a CSV of reference rows, or else a condition string for score_centiles()
  if (!is.null(opts$ref_data) && file.exists(opts$ref_data))
    opts$ref_data <- as.data.frame(data.table::fread(opts$ref_data))
  
  opts
}

#construct filename
model_filename <- function(pheno, split, total) {
  sprintf("%s_split%s_total%s_sharing_model.rds", pheno, split, as.character(total))
}

#locate one model: a file under model_dir if that is set, otherwise a URL into
#the repo at model_ref
model_path <- function(pheno, split, total, model_dir, model_ref) {
  f <- model_filename(pheno, split, total)
  if (!is.null(model_dir)) {
    file.path(model_dir, f)
  } else {
    sprintf("https://raw.githubusercontent.com/%s/%s/%s/%s",
            "BGDlab/sex_mod_braincharts", model_ref, "models_to_share", f)
  }
}

is_url <- function(p) grepl("^(https?|ftp)://", p)

read_model <- function(p) {
  if (!is_url(p)) return(readRDS(p))
  con <- gzcon(url(p, open = "rb"))
  on.exit(close(con), add = TRUE)
  readRDS(con)
}

#with a local model directory, say up if models are missing
check_model_dir <- function(model_dir, pheno_list, total, splits = c("A", "B")) {
  if (is.null(model_dir)) return(invisible(NULL))
  if (!dir.exists(model_dir))
    stop("model_dir does not exist: ", model_dir)

  want <- unlist(lapply(pheno_list, function(p)
    vapply(splits, function(s) model_path(p, s, total, model_dir, NULL), character(1))))
  gone <- want[!file.exists(want)]

  if (length(gone) == length(want))
    stop("no models for total=", as.character(total), " found in ", model_dir,
         "\n  expected files like ", basename(want[1]))
  if (length(gone) > 0)
    warning(length(gone), "/", length(want), " model file(s) not in ", model_dir,
            ":\n  ",
            paste(basename(utils::head(gone, 10)), collapse = "\n  "),
            if (length(gone) > 10) paste0("\n  ... and ", length(gone) - 10, " more"))

  invisible(gone)
}

#numeric training ranges, recovered from the model's pb() smoothers
training_ranges <- function(m) {
  rngs <- list()
  for (par in m$parameters) {
    sm <- colnames(m[[paste0(par, ".s")]])
    for (i in seq_along(sm)) {
      if (!grepl("^pb\\(", sm[i])) next
      v <- all.vars(str2lang(sm[i]))[1]
      #reaching into splinefun's environment is an internal, so fail soft
      r <- tryCatch(range(environment(gamlss::getSmo(m, par, which = i)$fun)$z$x),
                    error = function(e) NULL)
      if (length(r) != 2 || !is.numeric(r) || anyNA(r)) next
      rngs[[v]] <- if (is.null(rngs[[v]])) r else
        c(min(rngs[[v]][1], r[1]), max(rngs[[v]][2], r[2]))
    }
  }
  rngs
}

#find scorable rows for a given pheno
scorable_rows <- function(m, df, pheno, batch) {
  preds  <- unlist(list_predictors(m))
  needed <- c(pheno, preds)
  
  #check all cols available
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0) {
    warning(pheno, ": df is missing column(s) ", paste(missing_cols, collapse = ", "))
    return(integer(0))
  }
  
  #no NAs
  ok <- stats::complete.cases(as.data.frame(df)[, needed, drop = FALSE])
  
  #ignore batch variable
  check_preds <- setdiff(preds, batch)
  
  #check factor range - drop rows with unseen levels
  for (p in check_preds) {
    xlev <- unique(unlist(lapply(m$parameters, function(par) m[[paste0(par, ".xlevels")]][[p]])))
    if (length(xlev) == 0) next   #not a factor in the model
    unseen <- ok & !(as.character(df[[p]]) %in% xlev)
    if (any(unseen)) {
      warning(pheno, ": dropping ", sum(unseen), " row(s) with ", p,
              " level(s) not seen when the model was fit")
    }
    ok <- ok & !unseen
  }
  
  #check numeric range & warn if extrapolating
  rngs <- training_ranges(m)
  for (p in intersect(check_preds, names(rngs))) {
    v   <- suppressWarnings(as.numeric(df[[p]]))
    oob <- ok & !is.na(v) & (v < rngs[[p]][1] | v > rngs[[p]][2])
    if (any(oob)) {
      warning(pheno, ": ", sum(oob), " row(s) with ", p, " outside the fitted range [",
              signif(rngs[[p]][1], 6), ", ", signif(rngs[[p]][2], 6),
              "] - centiles are extrapolated")
    }
  }
  
  which(ok)
}

#stream one split's model, z-score the rows it can handle, then discard
score_split <- function(pheno, split, model_dir, model_ref, total, df, batch, ref_data, min_ref) {
  m <- tryCatch(
    read_model(model_path(pheno, split, total, model_dir, model_ref)), #local or github url
    error = function(e) {
      warning(pheno, " split ", split, ": could not read model (", conditionMessage(e), ")")
      NULL
    }
  )
  if (is.null(m)) return(NULL)
  on.exit({ rm(m); gc(verbose = FALSE) }, add = TRUE)
  
  idx <- scorable_rows(m, df, pheno, batch)
  if (length(idx) == 0) {
    warning(pheno, " split ", split, ": no scorable rows")
    return(NULL)
  }
  
  scores <- score_centiles(
    m,
    data        = as.data.frame(df)[idx, , drop = FALSE],
    fit_data    = NULL,
    standardize = TRUE,
    batch_term  = batch,
    ref_data = ref_data,
    min_ref = min_ref
  )
  
  out <- data.frame(.row_id = df$.row_id[idx], z = scores$std_score)
  names(out)[2] <- paste0("z_", split)
  out
}

#score a pheno against every split model and average the z-scores
score_pheno <- function(pheno, df, total, batch, splits = c("A", "B"), ref_data, min_ref,
                        model_dir, model_ref) {
  per_split <- Filter(Negate(is.null), lapply(splits, score_split,
                                              pheno = pheno,
                                              model_dir = model_dir,
                                              model_ref = model_ref,
                                              total = total, 
                                              df = df,
                                              batch = batch, 
                                              ref_data = ref_data, 
                                              min_ref = min_ref))
  if (length(per_split) == 0) {
    warning(pheno, ": no scores")
    return(NULL)
  }
  
  z_tbl  <- Reduce(function(a, b) dplyr::full_join(a, b, by = ".row_id"), per_split)
  z_cols <- setdiff(names(z_tbl), ".row_id")
  z_mat  <- as.matrix(z_tbl[, z_cols, drop = FALSE]) #pull zscores into matrix
  
  #average over whichever split models could score each subject
  n_splits <- rowSums(!is.na(z_mat))
  keep     <- n_splits > 0
  if (!any(keep)) {
    warning(pheno, ": no scorable subjects")
    return(NULL)
  }
  
  n_partial <- sum(n_splits[keep] < length(splits))
  if (n_partial > 0) {
    warning(pheno, ": ", n_partial, " subject(s) scored on fewer than ",
            length(splits), " split models")
  }
  
  z_mean <- rowMeans(z_mat[keep, , drop = FALSE], na.rm = TRUE)
  
  #reformat and convert mean back to centile space
  out <- data.frame(.row_id = z_tbl$.row_id[keep], z_mean, pnorm(z_mean))
  names(out) <- c(".row_id", paste0(pheno, c("_z", "_centile")))
  out
}
