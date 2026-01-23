#Fit GAMLSS models to select from on CV samples

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))
options(warn = 1)

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
base_mod <- readRDS(args[2])
save_path <- as.character(args[3])
total <- as.logical(args[4])
sim_df_uniform <- readRDS(args[5])

filename_no_ext <- sub("\\.[^.]*$", "", basename(args[2]))
filename <- sub("BestMod", "test", filename_no_ext)
file_full <- paste0(save_path, "/model_objs/", filename, "_full_mod.rds")

base_mod$call$data <- "df"
base_mod$call$family <- base_mod$family[[1]]

#FIT BASIC MODEL
#check if this pheno is already run, and if so, load
model <- NULL

if (file.exists(file_full)){
  print("loading pre-fit model")
  model <- tryCatch({readRDS(file_full)
  }, error = function(e){
    message(e$message, "- trying again")
    tryCatch({readRDS(file_full)
    }, error = function(e){
      message(e$message, "- refit")
      NULL
    })
  })
}
if (is.null(model)){
  model <- gamlss_lambda_rep(base_mod, null_mod="false", n.cyc=800)
}

#if model doesn't fit, fail
if (is.null(model)) {
  message("model fitting failed")
  stop("model fitting failed")
}

  #compare lambdas
  print("original lambdas:")
  print(base_mod$mu.lambda)
  print(base_mod$sigma.lambda)
  print(base_mod$nu.lambda)  

  print("refit lambdas:")
  print(model$mu.lambda)
  print(model$sigma.lambda)
  print(model$nu.lambda)
  
  mu_diff <- base_mod$mu.lambda - model$mu.lambda
  sig_diff <- base_mod$sigma.lambda - model$sigma.lambda
  stopifnot(c(mu_diff, sig_diff) < 0.01)
  if (!is.null(base_mod$nu.lambda)){
	nu_diff <- base_mod$nu.lambda - model$nu.lambda
	stopifnot(nu_diff < 0.01)}

print (paste("saving to", file_full))
saveRDS(model, file=file_full)

model$call$data <- "df"
model$call$family <- model$family[[1]]

#CENTILE FAN PLOT
pheno <- model$mu.terms[[2]] %>% as.character()

pred_list <- list_predictors(model)

#check if age is log-scaled
if ("logAge_days" %in% pred_list){
  age_var <- "logAge_days"
  print("simulate data for plotting")
  sim_df <- sim_data(df, age_var, factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")
  vars_of_interest <- c("sexMale_x_logAge", age_var, "sexMale")
} else {
  print("simulate data for plotting")
  age_var <- "age_days"
  sim_df <- sim_data(df, age_var, factor_var="sexMale", special_term = "sexMale_x_age = sexMale * age_days")
  vars_of_interest <- c("sexMale_x_age", age_var, "sexMale")
}

#define nuisance covars to residualize from points
resid_terms <- setdiff(pred_list, vars_of_interest)

birth <- log(280, base=10)
print("creating centile fan plot")
  fan_plot <- make_centile_fan(gamlssModel=model, 
                               df=df, 
                               x_var=age_var, 
                               color_var="sexMale",
                               get_peaks=TRUE, 
                               desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                               sim_data_list = sim_df,
                               remove_point_effect = resid_terms,
                               x_axis="log_lifespan_fetal")  +
    theme_linedraw() +
    labs(title=paste(pheno, "validation model"),
         x ="log Age (days)") +
  scale_color_discrete(name = "Sex", labels = c("Female", "Male")) +
    guides(fill=FALSE) +
    geom_vline(xintercept=birth)
    
    ggsave(file=paste0(save_path, "/centile_plots/", filename, ".png"), fan_plot)
    
  #WORM PLOT
    print("creating worm plot")
    wp <- wp.taki(xvar=df[[age_var]], resid=resid(model), n.inter=6)$plot +
      ggtitle(paste(pheno, "validation model"))
    ggsave(file=paste0(save_path, "/worm_plots/", filename, ".png"), wp)
    
  #COMPILE
    print("compiling stats")
    #centiles
    results_df <- cent_cdf(model, df, plot=FALSE, group="sexMale")
    
    #BIC & AIC
    summary_df <- data.frame(
      "AIC" = model$aic,
      "BIC" = model$sbc,
      "pheno" = pheno
    )

#SAVE CSVs
print("saving csvs")

#centiles
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", filename, "_centiles.csv"))

#BIC & AIC
fwrite(summary_df, file=paste0(save_path, "/model_sums/", filename, "_summary.csv"))

##################
#FIT NULL MODEL
print("fitting null model")

#check if null is already run, and if so, load
null_model <- NULL
file_null <- paste0(save_path, "/model_objs/", filename, "_null_mod.rds")

if (file.exists(file_null)){
  print("loading pre-fit model")
  null_model <- tryCatch({readRDS(file_null)
  }, error = function(e){
    message(e$message, "- trying again")
    tryCatch({readRDS(file_null)
    }, error = function(e){
      message(e$message, "- refit")
      NULL
    })
  })
}
if (is.null(null_model)){
  null_model <- gamlss_lambda_rep(base_mod, null_mod="true", n.cyc=800)
}

#if model doesn't fit, fail
if (is.null(null_model)) {
  message("null model fitting failed")
  stop("null model fitting failed")
}

print (paste("saving to", file_null))
saveRDS(null_model, file=file_null)

test_out <- LR.test(null_model, model, print=FALSE) #significance test
f2 <- cohens_f2_local(model, null_model) #effect size

#TEST
test_df <- data.frame(
  "chi" = test_out$chi,
  "df" = test_out$df,
  "p_val" = test_out$p.val,
  "fsq" = f2,
  "pheno" = pheno,
  "effect" = "sex_age"
)

##################
#FIT NULL MODEL WITH NO SEX-EFFECT FOR TBV-CORRECTED MODELS
if (total == TRUE) {
  print("fitting null model of all sex effects")
  null_model2 <- gamlss_lambda_rep(base_mod, null_mod="allSex", n.cyc=800)
  
  file_null2 <- paste0(save_path, "/model_objs/", filename, "_null2_mod.rds")
  print (paste("saving to", file_null2))
  saveRDS(null_model2, file=file_null2)
  
  test_out2 <- LR.test(null_model2, model, print=FALSE) #significance test
  f2_2 <- cohens_f2_local(model, null_model2) #effect size
  
  #TEST
  test_df2 <- data.frame(
    "chi" = test_out2$chi,
    "df" = test_out2$df,
    "p_val" = test_out2$p.val,
    "fsq" = f2_2,
    "pheno" = pheno,
    "effect" = "sex_all"
  )
  
  test_df <- rbind(test_df, test_df2)
}

fwrite(test_df, file=paste0(save_path, "/model_sums/", filename, "_LRtest.csv"))

####################
#GET SEX-DIFFERENCES
# Predict 50th cent & sigma values for each simulated level of factor_var

#initialize empty list(s)
centile_result_list <- list()
fname <- model$family[[1]]
qfun <- paste0("q", fname)
n_param <- length(model$parameters)

for (factor_level in names(sim_df)) {
  
  #make sure variable names are correct
  sub_df <- sim_df[[factor_level]]
  
  # Predict centiles
  print("predicting...")
  pred_df <- predictAll(model, newdata=sub_df, type="response", data=df)
  
  median_centile <- pred_centile(0.5, df = pred_df, q_func = qfun, n_param = n_param)
  centiles_df <- as.data.frame(median_centile)
  
  # check correct dim
  stopifnot(nrow(centiles_df) == nrow(pred_df))
  
  #add x_vals, predicted sigma name centiles for factor_var level and append to results list
  centiles_df$logAge_days <- sub_df$logAge_days
  
  # Sigma
  centiles_df$sigma <- pred_df$sigma
  sub_name <- paste0("pred_", factor_level)
  centile_result_list[[sub_name]] <- centiles_df
  
}

# zscore relative to male's initial y values - per taki
male_pheno <- df %>% 
  filter(sexMale==1) %>%
  pull(pheno)
mean_pheno <- mean(male_pheno)
std_dev_pheno <- sd(male_pheno)
z_score <- function(x){
  (x - mean_pheno) / std_dev_pheno
}

###for sigma, just use coefficient of variation (sigma over mean at each age)

#fun for derviative
get_deriv <- function(x, age) {
  c(NA, diff(x) / diff(age))
}

result_df <- bind_rows(centile_result_list, .id = "sexMale") %>%
  mutate(sex = if_else(sexMale=="pred_1", "Male", "Female")) %>%
  select(!sexMale) %>%
  tidyr:::pivot_wider(names_from=sex, values_from = c(median_centile, sigma)) %>%
  #z-score
  mutate(
    across(
      .cols = matches("centile"),
      .fns = z_score,
      .names = "{.col}_z"
    )) %>%
  #get coefficient of variation
  mutate(cv_Male = sigma_Male/median_centile_Male,
         cv_Female = sigma_Female/median_centile_Female) %>%
  #get_derivatives
  mutate(across(
    .cols = matches("centile|cv"),
    .fns = ~ get_deriv(.x, logAge_days),
    .names = "deriv_{.col}"
  )) %>%
  #get M-F differences
  mutate(centile_M_minus_F = median_centile_Male - median_centile_Female,
         centile_M_minus_F_z = median_centile_Male_z - median_centile_Female_z,
         cv_M_minus_F = cv_Male - cv_Female,
         deriv_M_minus_F = deriv_median_centile_Male - deriv_median_centile_Female,
         deriv_M_minus_F_z = deriv_median_centile_Male_z - deriv_median_centile_Female_z,
         deriv_CV_M_minus_F = deriv_cv_Male - deriv_cv_Female)

print("saving sex diffs")
fwrite(result_df, file=paste0(save_path, "/cent_csvs/", pheno, "_sexdiffs.csv"))

#### CALC ON SINGLE DATAFRAME ACROSS PHENOS ####
centile_result_list2 <- list()
min_age <- min(df$logAge_days)
max_age <- max(df$logAge_days)
for (factor_level in names(sim_df_uniform)) {
  
  #make sure variable names are correct
  sub_df <- sim_df_uniform[[factor_level]]
  sub_df <- sub_df %>%
    select(all_of(pred_list)) %>%
    filter(logAge_days >= min_age & logAge_days <= max_age)
  
  # Predict centiles
  print("predicting...")
  pred_df <- predictAll(model, newdata=sub_df, type="response", data=df)
  
  median_centile <- pred_centile(0.5, df = pred_df, q_func = qfun, n_param = n_param)
  centiles_df <- as.data.frame(median_centile)
  
  # check correct dim
  stopifnot(nrow(centiles_df) == nrow(pred_df))
  
  #add x_vals, predicted sigma name centiles for factor_var level and append to results list
  centiles_df$logAge_days <- sub_df$logAge_days
  
  # Sigma
  centiles_df$sigma <- pred_df$sigma
  sub_name <- paste0("pred_", factor_level)
  centile_result_list2[[sub_name]] <- centiles_df
  
}

result_df2 <- bind_rows(centile_result_list2, .id = "sexMale") %>%
  mutate(sex = if_else(sexMale=="pred_1", "Male", "Female")) %>%
  select(!sexMale) %>%
  tidyr:::pivot_wider(names_from=sex, values_from = c(median_centile, sigma)) %>%
  #z-score
  mutate(
    across(
      .cols = matches("centile|sigma"),
      .fns = z_score,
      .names = "{.col}_z"
    )) %>%
  #get_derivatives
  mutate(across(
    .cols = matches("centile"),
    .fns = ~ get_deriv(.x, logAge_days),
    .names = "deriv_{.col}"
  )) %>%
  #get M-F differences
  mutate(centile_M_minus_F = median_centile_Male - median_centile_Female,
         centile_M_minus_F_z = median_centile_Male_z - median_centile_Female_z,
         sigma_M_minus_F = sigma_Male - sigma_Female,
         sigma_M_minus_F_z = sigma_Male - sigma_Female_z,
         deriv_M_minus_F = deriv_median_centile_Male - deriv_median_centile_Female,
         deriv_M_minus_F_z = deriv_median_centile_Male_z - deriv_median_centile_Female_z)

print("saving sex diffs")
fwrite(result_df2, file=paste0(save_path, "/cent_csvs/", pheno, "_uniform_sexdiffs.csv"))

print("SUCCESS")
