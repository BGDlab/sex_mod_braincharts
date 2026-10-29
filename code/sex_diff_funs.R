#Functions for estimating sex differences from fit GAMLSS models
#requires dplyr, tidyr, gamlss (incl. gamlss.dist), gamlssTools

#fun for derviative
get_deriv <- function(x, age) {
  c(NA, diff(x) / diff(age))
}

# Function to process simulated datasets and calculate sex differences
# Predict 50th cent, sigma & centile-based SD values for each simulated level of sexMale
# M-F differences are z-scored against the male centile-based SD at each age (_agez)
#   sim_data_list: list of simulated dfs, one per level of sexMale (e.g. from sim_data())
#   model: fit gamlss model
#   df: data model was fit on
#   age_var: name of age column in model ("age_days" or "logAge_days")
#   pred_list: predictors in model (from list_predictors())
#   pheno: name of phenotype column in df (no longer used; kept so existing calls still work)
process_sex_diffs <- function(sim_data_list, model, df, age_var, pred_list, pheno) {

  fname <- model$family[[1]]
  qfun <- paste0("q", fname)
  n_param <- length(model$parameters)

  # Initialize empty list
  centile_result_list <- list()

  # Get age range of fit data to filter sim grid
  min_age <- min(df[[age_var]])
  max_age <- max(df[[age_var]])

  # Process each factor level
  for (factor_level in names(sim_data_list)) {
    sub_df <- sim_data_list[[factor_level]] %>%
      select(all_of(pred_list)) %>%
      filter(.data[[age_var]] >= min_age & .data[[age_var]] <= max_age)

    # Predict centiles
    print("predicting...")
    pred_df <- predictAll(model, newdata=sub_df, type="response", data=df)

    median_centile <- pred_centile(0.5, df = pred_df, q_func = qfun, n_param = n_param)
    centiles_df <- as.data.frame(median_centile)

    # check correct dim
    stopifnot(nrow(centiles_df) == nrow(pred_df))

    #add x_vals, predicted sigma name centiles for factor_var level and append to results list
    centiles_df[[age_var]] <- sub_df[[age_var]]

    # Sigma
    centiles_df$sigma <- pred_df$sigma

    # Centile-based SD at each age, from the fitted distribution (IQR / 1.349)
    centiles_df$sd <- centile_sd(fname, pred_df$mu, pred_df$sigma, pred_df$nu, pred_df$tau)
    sub_name <- paste0("pred_", factor_level)
    centile_result_list[[sub_name]] <- centiles_df
  }

  # Process results
  result_df <- bind_rows(centile_result_list, .id = "sexMale") %>%
    mutate(sex = if_else(sexMale=="pred_1", "Male", "Female")) %>%
    select(!sexMale) %>%
    tidyr:::pivot_wider(names_from=sex, values_from = c(median_centile, sigma, sd))

  result_df <- result_df %>%
    #get ratios
    mutate(logcv_M_div_F = log(sigma_Male/sigma_Female),
           log_M_div_F = log(median_centile_Male/median_centile_Female)) %>%
    #get derivs
    mutate(across(
      .cols = matches("centile|_div_F"),
      .fns = ~ get_deriv(.x, .data[[age_var]]),
      .names = "deriv_{.col}"
    )) %>%
  #get M-F differences
    mutate(centile_M_minus_F = median_centile_Male - median_centile_Female,
           deriv_M_minus_F = deriv_median_centile_Male - deriv_median_centile_Female,
           #z score against the male centile-based SD at each age
           centile_M_minus_F_agez = if_else(is.na(sd_Male), NA_real_, centile_M_minus_F / sd_Male),
           #derivative of the standardized difference itself: with an age-varying SD, this isn't
           #the difference of the sexes' derivatives divided by the SD
           deriv_M_minus_F_agez = get_deriv(centile_M_minus_F_agez, .data[[age_var]]))

  return(result_df)
}

# Mean and SD of a gamlss.dist distribution, from its density:
#   E(Y) = integral of y f(y), Var(Y) = integral of [y - E(Y)]^2 f(y), SD(Y) = sqrt(Var(Y))
#   family: gamlss family name, e.g. "BCTo", "GG", "NO", "PO"
#   params: named list of parameter values, e.g. list(mu = 100, sigma = 0.1, nu = 0.5, tau = 10)
# Returns NA (with a warning) when the mean or variance is not finite, e.g. heavy-tailed BCTo
dist_moments <- function(family, params) {
  fam <- gamlss.dist::as.gamlss.family(family)
  dfun <- get(paste0("d", family), mode = "function")
  qfun <- get(paste0("q", family), mode = "function")
  f <- function(y) do.call(dfun, c(list(y), params))
  q <- function(p) do.call(qfun, c(list(p), params))

  if (fam$type == "Discrete") {
    #sum over the support 0, 1, 2, ... up to where the remaining probability is negligible
    y <- seq(0, q(1 - 1e-12))
    m1 <- sum(y * f(y))
    v <- sum((y - m1)^2 * f(y))
  } else if (fam$type == "Continuous") {
    #support: whole real line, positive values, or (0, 1)
    lower <- if (fam$y.valid(-1)) -Inf else 0
    upper <- if (fam$y.valid(2)) Inf else 1
    #split at central quantiles so integrate() finds the mass wherever it sits
    breaks <- c(lower, q(c(0.001, 0.5, 0.999)), upper)
    piece_int <- function(g) {
      sum(sapply(seq_len(length(breaks) - 1),
                 function(i) integrate(g, breaks[i], breaks[i + 1])$value))
    }
    #integrate() fails when the mean or variance is infinite (heavy tails), so return NA
    out <- tryCatch({
      m1 <- piece_int(function(y) y * f(y))
      v <- piece_int(function(y) (y - m1)^2 * f(y))
      c(m1, v)
    }, error = function(e) {
      warning("mean/variance not finite for ", family, "(",
              paste(names(params), signif(unlist(params), 3), sep = "=", collapse = ", "),
              "): ", conditionMessage(e), call. = FALSE)
      c(NA, NA)
    })
    m1 <- out[1]
    v <- out[2]
  } else {
    stop("mixed distributions (e.g. zero-adjusted) are not supported")
  }
  c(mean = m1, sd = sqrt(v))
}

# Vectorized SD, e.g. one value per age from predictAll() output
dist_sd <- function(family, mu, sigma, nu = NULL, tau = NULL) {
  params <- Filter(Negate(is.null), list(mu = mu, sigma = sigma, nu = nu, tau = tau)) #drop unused moments
  do.call(mapply, c(list(FUN = function(...) dist_moments(family, list(...))[["sd"]]), params))
}

# Centile-based approximation of the SD, e.g. one value per age from predictAll() output
# Half the distance between the p and 1-p centiles, scaled so it equals the SD for a normal
# distribution: p = 0.75 gives IQR / 1.349 (the semi-interquartile range, Rigby et al. 2019
# sec 2.3), p = pnorm(1) gives half the 16th-84th centile range
# Unlike dist_sd(), it always exists, even when heavy tails make the variance infinite
centile_sd <- function(family, mu, sigma, nu = NULL, tau = NULL, p = 0.75) {
  qfun <- get(paste0("q", family), mode = "function")
  params <- Filter(Negate(is.null), list(mu = mu, sigma = sigma, nu = nu, tau = tau)) #drop unused moments
  (do.call(qfun, c(list(p), params)) - do.call(qfun, c(list(1 - p), params))) / (2 * qnorm(p))
}
