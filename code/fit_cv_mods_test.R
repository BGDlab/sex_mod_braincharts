#Fit GAMLSS models to select from on CV samples

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

source("./code/gamlss_fit_funs.R")

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
pheno <- as.character(args[2])
base_mod <- readRDS(args[3])
save_path <- as.character(args[4])
weight_pts <- as.logical(args[5])
log_scale <- as.logical(args[6])

#drop extra variables
df <- df %>%
  dplyr::select(all_of(c(pheno, fs, "logAge_days", "sexMale", "study_site", "sexMale_x_logAge", "age_days"))) %>%
  na.omit()

#inverse-weight by age w/in sex (written w help from gpt)
if (weight_pts == TRUE){
  n_bins <- 50
  df$age_bin <- cut(df$age_days, breaks = n_bins, include.lowest = TRUE)
  
  
  df <- df %>%
    mutate(sex = as.factor(sexMale)) %>%
    group_by(sex, age_bin) %>%
    mutate(bin_count = n()) %>%
    group_by(sex) %>%
    mutate(
      observed_prob = bin_count / sum(bin_count),
      uniform_prob = 1 / n_bins,
      raw_weight = uniform_prob / observed_prob,
      weight = raw_weight * (n() / sum(raw_weight))
    ) %>%
    ungroup() %>%
    select(-bin_count, -observed_prob, -uniform_prob, - raw_weight)
}

#log-scale pheno if necessary
if (log_scale == TRUE){
  pheno_sym <- sym(pheno)
  
  df <- df %>%
    mutate(!!pheno_sym := ifelse(!!pheno_sym==0, 1, !!pheno_sym)) %>% #replace 0 with 1
    mutate(!!pheno_sym := log(!!pheno_sym, base=10)) #transform
}

#sim data ONCE for centile fan plotting
print("simulate data for plotting")
sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")

#FIT BASIC MODEL
  model <- gamlss_3lambda_rep(base_mod, null_mod=FALSE)
    

  #if model isn't fit, skip to next loop
  if (is.null(model)) {
    message("model fitting failed")
    stop()
  }
   
  saveRDS(model, file=paste0(save_path, "/model_objs/", pheno, "_validate_full_mod.rds"))
    
  stopifnot(base_mod$mu.lambda == model$mu.lambda)
  stopifnot(base_mod$sigma.lambda == model$sigma.lambda)
  
  #CENTILE FAN PLOT
    print("creating centile fan plot")
    fan_plot <- make_centile_fan(gamlssModel=model, df=df, x_var="logAge_days", color_var="sexMale",
                               get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                               sim_data_list = sim_df) +
    ggtitle(paste(pheno, "validation model"))
    
    ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_validation.png"), fan_plot)
    
  #WORM PLOT
    print("creating worm plot")
    wp <- wp.taki(xvar=df$logAge_days, resid=resid(model), n.inter=8) +
      ggtitle(paste(pheno, "validation model"))
    ggsave(file=paste0(save_path, "/worm_plots/", pheno, "_validation.png"), wp)
    
  #COMPILE
    print("compiling stats")
    #centiles
    results_df <- cent_cdf(model, df, "sexMale")
    results_df$lambda <- l.name
    
    #BIC & AIC
    summary_df <- data.frame(
      "AIC" = model$aic,
      "BIC" = model$sbc,
      "lambda" = l.name,
      "pheno" = pheno,
      "fs_moment" = fs_include
    )

#SAVE CSVs
print("saving csvs")

#centiles
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", pheno, "_validation_results.csv"))

#BIC & AIC
fwrite(summary_df, file=paste0(save_path, "/model_sums/", pheno, "_validation_summary.csv"))

##################
#FIT NULL MODEL
null_model <- gamlss_3lambda_rep(base_mod, null_mod=TRUE)
test_out <- LR.test(null_model, model, print=FALSE)

#TEST
test_df <- data.frame(
  "AIC" = test_out$chi,
  "BIC" = test_out$sbc,
  "p_val" = test_out$p.val,
  "pheno" = pheno
)
fwrite(test_df, file=paste0(save_path, "/model_sums/", pheno, "_test.csv"))

