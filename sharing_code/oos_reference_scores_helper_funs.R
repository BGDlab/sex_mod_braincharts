#helper functions for getting out-of-sample reference scores

model_url <- function(pheno, split, total) {
  sprintf(
    "https://raw.githubusercontent.com/%s/%s/%s/%s_split%s_total%s_sharing_model.rds",
    "BGDlab/sex_mod_braincharts", model_ref, "models_to_share", pheno, split, as.character(total)
  )
}

read_rds_url <- function(u) {
  con <- gzcon(url(u, open = "rb"))
  on.exit(close(con), add = TRUE)
  readRDS(con)
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
score_split <- function(pheno, split, total, df, batch) {
  m <- tryCatch(
    read_rds_url(model_url(pheno, split, total)),
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
    batch_term  = batch
  )
  
  out <- data.frame(.row_id = df$.row_id[idx], z = scores$std_score)
  names(out)[2] <- paste0("z_", split)
  out
}

#score a pheno against every split model and average the z-scores
score_pheno <- function(pheno, df, total, batch, splits = c("A", "B")) {
  per_split <- Filter(Negate(is.null), lapply(splits, score_split,
                                              pheno = pheno, total = total, df = df,
                                              batch = batch))
  if (length(per_split) == 0) {
    warning(pheno, ": no scores")
    return(NULL)
  }
  
  z_tbl  <- Reduce(function(a, b) dplyr::full_join(a, b, by = ".row_id"), per_split)
  z_cols <- setdiff(names(z_tbl), ".row_id")
  z_mat  <- as.matrix(z_tbl[, z_cols, drop = FALSE])
  
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
  
  out <- data.frame(.row_id = z_tbl$.row_id[keep], z_mean, pnorm(z_mean))
  names(out) <- c(".row_id", paste0(pheno, c("_z", "_centile")))
  out
}
