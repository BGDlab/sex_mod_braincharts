#Functions for estimating sex differences from fit GAMLSS models
#requires dplyr, tidyr, gamlss, gamlssTools

#fun for derviative
get_deriv <- function(x, age) {
  c(NA, diff(x) / diff(age))
}

# Function to process simulated datasets and calculate sex differences
# Predict 50th cent & sigma values for each simulated level of sexMale
#   sim_data_list: list of simulated dfs, one per level of sexMale (e.g. from sim_data())
#   model: fit gamlss model
#   df: data model was fit on
#   age_var: name of age column in model ("age_days" or "logAge_days")
#   pred_list: predictors in model (from list_predictors())
#   pheno: name of phenotype column in df
process_sex_diffs <- function(sim_data_list, model, df, age_var, pred_list, pheno) {

  # zscore relative to male's initial y values - per taki
  male_pheno <- df %>%
    filter(sexMale==1) %>%
    pull(all_of(pheno))
  mean_pheno <- mean(male_pheno)
  std_dev_pheno <- sd(male_pheno)
  z_score <- function(x){
    (x - mean_pheno) / std_dev_pheno
  }

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
    sub_name <- paste0("pred_", factor_level)
    centile_result_list[[sub_name]] <- centiles_df
  }

  # Process results
  result_df <- bind_rows(centile_result_list, .id = "sexMale") %>%
    mutate(sex = if_else(sexMale=="pred_1", "Male", "Female")) %>%
    select(!sexMale) %>%
    tidyr:::pivot_wider(names_from=sex, values_from = c(median_centile, sigma))

  result_df <- result_df %>%
    #z score
    mutate(
      across(
        .cols = matches("centile"),
        .fns = z_score,
        .names = "{.col}_z"
      )) %>%
    #get CV ratios
    mutate(cv_M_div_F = sigma_Male/sigma_Female,
           logcv_M_div_F = log(sigma_Male/sigma_Female)) %>%
    #get derivs
    mutate(across(
      .cols = matches("centile|cv"),
      .fns = ~ get_deriv(.x, .data[[age_var]]),
      .names = "deriv_{.col}"
    )) %>%
  #get M-F differences
    mutate(centile_M_minus_F = median_centile_Male - median_centile_Female,
           centile_M_minus_F_z = median_centile_Male_z - median_centile_Female_z,
           deriv_M_minus_F = deriv_median_centile_Male - deriv_median_centile_Female,
           deriv_M_minus_F_z = deriv_median_centile_Male_z - deriv_median_centile_Female_z)

  return(result_df)
}
