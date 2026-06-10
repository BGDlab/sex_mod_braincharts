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
full_df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
fs <- as.character(args[2]) #freesurfer version
base_mod <- readRDS(args[3])
save_path <- as.character(args[4])
total <- as.character(args[5])

filename_no_ext <- sub("\\.[^.]*$", "", basename(args[3]))
filename <- sub("BestMod", "train", filename_no_ext)
file_full <- paste0(save_path, "/model_objs/", filename, "_age2plus.rds")

#check if this pheno is already run, and if so, end
if (file.exists(file_full)){
stop("Already tested, skipping pheno")
}

##### READ INFO #####
base_mod$call$data <- "df"
pheno <- base_mod$mu.terms[[2]] %>% as.character()
pred_list <- list_predictors(base_mod)

#check if age is log-scaled
if ("logAge_days" %in% pred_list){
  age_var <- "logAge_days"
  sex_age_var <- "sexMale_x_logAge"
} else {
  age_var <- "age_days"
  sex_age_var <- "sexMale_x_age"
}
#define nuisance covars to residualize from points
vars_of_interest <- c(age_var, sex_age_var, "sexMale")
resid_terms <- setdiff(pred_list, vars_of_interest)

##### PREP DATAFRAME #####
#drop extra variables
if (total == "FALSE"){
  df <- full_df %>%
    dplyr::select(all_of(c(pred_list, pheno))) %>%
    na.omit() %>%
    trunc_coverage(age_var, max_loops=100) #drop points at ends if too sparse
} else {
  df <- full_df %>%
    dplyr::select(all_of(c(pred_list, pheno))) %>%
    na.omit() %>%
    trunc_coverage(c(total, age_var), max_loops=100) #drop points at ends if too sparse
}

##### FIT TRAINING MODEL #####
model <- gamlss_lambda_rep(base_mod, 
                           null_mod="false",
                           keep_lambdas=FALSE,
                           n.cyc=1000)

#if model isn't fit, skip to next loop
if (is.null(model)) {
  message("model fitting failed")
  stop("model fitting failed")
}

print (paste("saving to", file_full))
saveRDS(model, file=file_full)

model$call$data <- "df"
model$call$family <- model$family[[1]]

#CENTILE FAN PLOT
print("creating centile fan plot")
#check if age is log-scaled
if ("logAge_days" %in% pred_list){
  print("simulate data for plotting")
  sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")
} else {
  print("simulate data for plotting")
  sim_df <- sim_data(df, "age_days", factor_var="sexMale", special_term = "sexMale_x_age = sexMale * age_days")
}

  fan_plot <- make_centile_fan(gamlssModel=model, 
                               df=df, 
                               x_var=age_var, 
                               color_var="sexMale",
                               get_peaks=TRUE, 
                               desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                               sim_data_list = sim_df,
                               remove_point_effect = resid_terms)  +
    labs(title=paste(pheno, "model, age 2+ yrs"),
         x ="log Age (days)",
         color = "Sex=Male", fill="Sex=Male")
    
    ggsave(file=paste0(save_path, "/centile_plots/", filename, ".png"), fan_plot)
    
  #WORM PLOT
    print("creating worm plot")
    wp <- wp.taki(xvar=df$logAge_days, resid=resid(model), n.inter=6)$plot +
      ggtitle(paste(pheno, "model, age 2+ yrs"))
    ggsave(file=paste0(save_path, "/worm_plots/", filename, ".png"), wp)
    
  #COMPILE
    print("compiling stats")
    #centiles
    #print(names(df))
    #print(list_predictors(model))
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

print("SUCCESS")
