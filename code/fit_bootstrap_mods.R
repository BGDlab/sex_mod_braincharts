#bootstrap brainchart model & save predicted trajectories
#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))

args <- commandArgs(trailingOnly = TRUE)
print(args)
n <- args[1]
pheno <- as.character(args[2])
og_df <- fread(args[3], stringsAsFactors = TRUE, na.strings = "") #path to csv
base_mod <- readRDS(args[4])
save_path <- as.character(args[5])

#### BOOTSTRAP DF ####
df <- og_df %>%
  mutate(study = as.factor(sub("_.*", "", study_site)), 
         sex = as.factor(sexMale)) %>%
  slice_sample(prop=1, by=c(sex, study), replace=TRUE)

#### FIT MODEL ####
model <- gamlss_lambda_rep(base_mod,
                           null_mod=FALSE,
                           keep_lambdas=FALSE,
                           start.from=NULL,
                           weight=FALSE)

#### PREDICT OUTPUTS 
print("predicting...")

#check if age is log-scaled
sim_df <- sim_data(og_df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")

#initialize empty list(s)
centile_result_list <- list()
fname <- model$family[[1]]
qfun <- paste0("q", fname)
n_param <- length(model$parameters)


# Predict 50th cent & sigma values for each simulated level of factor_var
for (factor_level in names(sim_df)) {
  
  #make sure variable names are correct
  sub_df <- sim_df[[factor_level]]
  
  # Predict centiles
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

deriv_df <- data.frame(deriv_Male = male_deriv,
                       deriv_Female = female_deriv,
                       logAge_days = x_mid)

final_df <- full_join(result_df, deriv_df) %>%
  mutate(centile_M_minus_F = median_centile_Male - median_centile_Female,
         sigma_M_minus_F = sigma_Male - sigma_Female,
         deriv_M_minus_F = deriv_Male - deriv_Female)

final_df$pheno <- pheno
final_df$boot <- n

print("saving sex diffs")

#### SAVE ####
fwrite(final_df, file=paste0(save_path, "/cent_csvs/", pheno, "_boot", n, "_out.csv"))
