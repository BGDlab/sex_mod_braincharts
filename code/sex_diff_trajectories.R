#calculate male-female trajectory differences, using code from fit_brainchart_mods.R
set.seed(555566)
#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))

#read in model & dataframe
#test models are saved in the split that they're tested on (i.e. cv_sample_A_test has models that were trained on B)
args <- commandArgs(trailingOnly = TRUE)
print(args)
model <- readRDS(args[1])
df <- fread(args[2], stringsAsFactors = TRUE, na.strings = "") #path to csv
save_path <- as.character(args[3])

#initialize empty list(s)
centile_result_list <- list()

#setup
fname <- model$family[[1]]
qfun <- paste0("q", fname)
n_param <- length(model$parameters)
model$call$data <- "df"
model$call$family <- fname
pred_list <- list_predictors(model)

if ("logAge_days" %in% pred_list){
  age_var <- "logAge_days"
  print("simulate data for plotting")
  sim_df <- sim_data(df, age_var, factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")
} else {
  print("simulate data for plotting")
  age_var <- "age_days"
  sim_df <- sim_data(df, age_var, factor_var="sexMale", special_term = "sexMale_x_age = sexMale * age_days")
}

# Predict 50th cent & sigma values for each simulated level of factor_var
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

result_df <- bind_rows(centile_result_list, .id = "sexMale") %>%
  mutate(sex = if_else(sexMale=="pred_1", "Male", "Female")) %>%
  select(!sexMale) %>%
  tidyr:::pivot_wider(names_from=sex, values_from = c(median_centile, sigma))

#Derivative diffs
dy_male <- diff(result_df$median_centile_Male)
dy_female <- diff(result_df$median_centile_Female)
dx <- diff(result_df$logAge_days)

male_deriv <- dy_male/dx
female_deriv <- dy_female/dx

x_mid <- zoo::rollmean(result_df$logAge_days, 2)

deriv_df <- data.frame("deriv_Male" = male_deriv,
                       "deriv_Female" = female_deriv,
                       "logAge_days" = x_mid)

final_df <- full_join(result_df, deriv_df) %>%
  mutate(centile_M_minus_F = median_centile_Male - median_centile_Female,
         sigma_M_minus_F = sigma_Male - sigma_Female,
         deriv_M_minus_F = deriv_Male - deriv_Female)

print("saving sex diffs")
fwrite(final_df, file=paste0(save_path, "/cent_csvs/", pheno, "_sexdiffs.csv"))