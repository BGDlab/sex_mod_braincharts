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
  } else {
    mod_to_fit <- mod_B
  }

} else if (length(term_list_A) < length(term_list_B)) {
  # use mod A
  print("using model A")
  mod_to_fit <- mod_A
  
} else if (length(term_list_A) > length(term_list_B)){
  # use mod B
  print("using model B")
  mod_to_fit <- mod_B
  
} else {
  stop("can't find simpler model")
}

######### FIT MODEL #########

model <- gamlss_lambda_rep(mod_to_fit,
                   null_mod=FALSE,
                   keep_lambdas=FALSE,
                   start.from=NULL,
                   weight=FALSE)

#if model isn't fit, skip to next loop
if (is.null(model)) {
  message("model fitting failed")
  stop()
}

print(paste("saving", pheno, "model"))
file_full <- paste0(save_path, "/model_objs/", pheno, "_brainchart.rds")
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
