#Find simpler model between samples A and B
set.seed(555566)
#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))

#load pheno, lbcc df, cv_A and cv_B mods
args <- commandArgs(trailingOnly = TRUE)
print(args)
pheno <- as.character(args[1])
df <- fread(args[2], stringsAsFactors = TRUE, na.strings = "") #path to csv
mod_A <- readRDS(args[3])
mod_B <- readRDS(args[4])
save_path <- as.character(args[5])

file_full <- paste0(save_path, "/model_objs/", pheno, "_brainchart.rds")
#check if this pheno is already run, and if so, end
if (file.exists(file_full)){
  stop("Brainchart already fit, skipping pheno")
}

#load both A and B models
mod_A$call$data <- "df"
mod_B$call$data <- "df"

cv_mod_list <- list("A" = list(mod_A),
                    "B" = list(mod_B))

#list predictors in each model/moment
term_list_A <- c(list_predictors(mod_A, moment="mu"), list_predictors(mod_A, moment="sigma"), list_predictors(mod_A, moment="nu"))
term_list_B <- c(list_predictors(mod_B, moment="mu"), list_predictors(mod_B, moment="sigma"), list_predictors(mod_B, moment="nu"))

if (length(term_list_A) == length(term_list_B)) {
  #if lengths are equal, check if identical models
  if (all(rm_lambdas(mod_A$mu.formula) == rm_lambdas(mod_B$mu.formula)) &
      all(rm_lambdas(mod_A$sigma.formula) == rm_lambdas(mod_B$sigma.formula)) &
      all(rm_lambdas(mod_A$nu.formula) == rm_lambdas(mod_B$nu.formula))
      ){
    print("models from samples A and B are identical")
  } else {
    print("models from A and B are NOT identical")
  }
  
  #pick random mod
  random_mod <- sample(c("A", "B"), 1)
  print(paste("randomly selecting model", random_mod))
  if (random_mod == "A"){
    mod_to_fit <- mod_A
    alt_mod <- mod_B
  } else {
    mod_to_fit <- mod_B
    alt_mod <- mod_A
  }

} else if (length(term_list_A) < length(term_list_B)) {
  # use mod A
  print("using model A")
  mod_to_fit <- mod_A
  alt_mod <- mod_B
  
} else if (length(term_list_A) > length(term_list_B)){
  # use mod B
  print("using model B")
  mod_to_fit <- mod_B
  alt_mod <- mod_A
  
} else {
  stop("can't find simpler model")
}

#clean up
rm(mod_A)
rm(mod_B)

######### FIT MODEL #########

model <- gamlss_lambda_rep(mod_to_fit,
                   null_mod="false",
                   keep_lambdas=FALSE,
                   start.from=NULL,
                   weight=NULL)

#if model isn't fit, try other model
if (is.null(model)) {
  message("model fitting failed, trying other CV mod")
  model <- gamlss_lambda_rep(alt_mod,
                             null_mod="false",
                             keep_lambdas=FALSE,
                             start.from=NULL,
                             weight=NULL)
}

#if fails again, stop
if (is.null(model)) {
  message("model fitting failed")
  stop()
}


print(paste("saving", pheno, "model"))
saveRDS(model, file=file_full)

model$call$data <- "df"
model$call$family <- model$family[[1]]

######### PLOTS #########
birth <- log(280, base=10)
pred_list <- list_predictors(model)
#see if fs_version included as covariate
fs <- pred_list[grep("^fs_version", pred_list)]
if (length(fs) == 1){
  resid_terms <- c(fs, "study_site")
} else {
  resid_terms <- "study_site"
}

#check if age is log-scaled
if ("logAge_days" %in% pred_list){
  age_var <- "logAge_days"
  print("simulate data for plotting")
  sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")
} else {
  print("simulate data for plotting")
  sim_df <- sim_data(df, "age_days", factor_var="sexMale", special_term = "sexMale_x_age = sexMale * age_days")
  age_var <- "age_days"
}


#CENTILE FAN PLOT
print("creating centile fan plot")
fan_plot <- make_centile_fan(gamlssModel=model,
                             df=df, 
                             x_var=age_var, 
                             color_var="sexMale",
                             get_peaks=TRUE, 
                             desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                             sim_data_list = sim_df,
                             remove_point_effect = resid_terms,
                             x_axis="log_lifespan_fetal",
                             y_scale=unscale)  +
  theme_linedraw() +
  scale_color_discrete(name = "Sex", labels = c("Female", "Male")) +
  guides(fill=FALSE) +
  geom_vline(xintercept=birth)
  labs(title=paste(pheno, "brainchart"),
       x ="Age (log days)")

ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_centilefan.png"), fan_plot)

#WORM PLOT
print("creating worm plot")
wp <- wp.taki(xvar=df$logAge_days, resid=resid(model), n.inter=10) 
plt <- wp$plot
ggsave(file=paste0(save_path, "/worm_plots/", pheno, "_wp.png"), plt)
print(wp$outliers)
fwrite(wp$outliers, file=paste0(save_path, "/worm_plots/", pheno, "_outliers.csv"))

#COMPILE
print("compiling stats")
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
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", pheno, "_centiles.csv"))

#BIC & AIC
fwrite(summary_df, file=paste0(save_path, "/model_sums/", pheno, "_summary.csv"))

### MALE - FEMALE DIFFS
print("calculating sex diffs")

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
