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
l.name <- as.character(args[3])
fs <- as.character(args[4])
fs_include <- as.character(args[5])
save_path <- as.character(args[6])
log_scale <- as.logical(args[7])

#drop extra variables
df <- df %>%
  dplyr::select(all_of(c(pheno, fs, "logAge_days", "sexMale", "study_site", "sexMale_x_logAge"))) %>%
  na.omit()

#log-scale pheno if necessary
if (log_scale == TRUE){
  pheno_sym <- sym(pheno)
  
  df <- df %>%
    mutate(!!pheno_sym := ifelse(!!pheno_sym==0, 1, !!pheno_sym)) %>% #replace 0 with 1
    mutate(!!pheno_sym := log(!!pheno_sym, base=10)) #transform
}

#define lambdas to be tested
if (l.name == "NULL"){
  l <- NULL
} else {
  l <- as.numeric(l.name)
}

results_df <- data.frame()
summary_df <- data.frame()

#sim data ONCE for centile fan plotting
print("simulate data for plotting")
sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")

#FIT MODEL
print(paste("fitting model with lambda =", l, "and fs in", fs_include))
  
#FIT BASIC MODEL
  model <- gamlss_3lambda(pheno, lambda=l, fs_ver=fs, fs_moment=fs_include, fam="GG")
    

  #if model isn't fit, skip to next loop
  if (is.null(model)) {
    message("model fitting failed")
    stop()
  }
   
  saveRDS(model, file=paste0(save_path, "/model_objs/", pheno, "_lambda", l.name, "_", fs_include, "_mod.rds"))
    
  #CENTILE FAN PLOT
    print("creating centile fan plot")
    fan_plot <- make_centile_fan(gamlssModel=model, df=df, x_var="logAge_days", color_var="sexMale",
                               get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                               sim_data_list = sim_df) +
    ggtitle(paste(pheno, "\nsmoothed w/ lamda=", l.name))
    
    ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_lambda", l.name, "_", fs_include, ".png"), fan_plot)
    
  #WORM PLOT
    print("creating worm plot")
    wp <- wp.taki(xvar=df$logAge_days, resid=resid(model), n.inter=8) +
      ggtitle(paste(pheno, "\nsmoothed w/ lambda=", l.name))
    ggsave(file=paste0(save_path, "/worm_plots/", pheno, "_lambda", l.name, "_", fs_include, ".png"), wp)
    
  #COMPILE
    print("compiling stats")
    #centiles
    results_df <- cent_cdf(model, df, "sexMale")
    results_df$lambda <- l.name
    
    #results_df <- rbind(results_df, sub_df)
    
    #BIC & AIC
    summary_df <- data.frame(
      "AIC" = model$aic,
      "BIC" = model$sbc,
      "lambda" = l.name,
      "pheno" = pheno
    )
    
    #summary_df <- rbind(summary_df, sum_df)
    
    #z-score normality - NEED TO FIGURE THIS OUT
    #Q.stats(model, xvar=logAge_days, n.inter=5, plot=FALSE)


#SAVE CSVs
print("saving csvs")

#centiles
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", pheno, "_lambda", l.name, "_", fs_include, "_results.csv"))

#BIC & AIC
fwrite(summary_df, file=paste0(save_path, "/model_sums/", pheno, "_lambda", l.name, "_", fs_include, "_summary.csv"))

#z-score normality

###################
#print(paste(sum(unique(results_df$lambda)), "of", loop_count, pheno, fs_include, "models successful"))
