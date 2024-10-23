#sub-divide IDPs into the largest possible complete chunks for ComBat harmonization

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

source("gamlss_fit_funs.R")

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
pheno <- as.character(args[2])
knot_lists <- readRDS(as.character(args[3]))
save_path <- as.character(args[4])
log_scale <- as.logical(args[5])
fs_covary <- as.logical(args[6])

#drop extra variables
df <- df %>%
  dplyr::select(all_of(c(pheno, "logAge_days", "sexMale", "fs_version", "study_site", "sexMale_x_logAge"))) %>%
  na.omit()

#log-scale pheno if necessary
if (log_scale == TRUE){
  df <- df %>%
    mutate(pheno = ifelse(pheno==0, 1, pheno)) %>% #replace 0 with 1
  mutate(across(c(pheno), \(x) log(x, base=10))) #transform
}

#define degrees of freedom to be tested
degree_list <- seq(3,6)

results_df <- data.frame()

#sim data ONCE for centile fan plotting
print("simulate data for plotting")
sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")

loop_count <- 0

#FIT BASE MODEL
for (degree in degree_list){
  
  knot_index <- paste0("df", degree)
  knots_list <- paste(knot_lists[[knot_index]], collapse=", ")

  for (sigma_degree in c(ceiling(degree/2), degree)){
   if (sigma_degree <2){
     next
   }
    s_knot_index <- paste0("df", sigma_degree)
    s_knots_list <- paste(knot_lists[[s_knot_index]], collapse=", ")
    
    print(paste("fitting model with df = ", degree, "in mu and df =", sigma_degree, "in sigma"))
    
    #fit model with or without fs_version term
    if (fs_covary==TRUE){
      model <- gamlss_mod_nofs(pheno, knots=knots_list, sigma_knots=s_knots_list)
    } else {
      model <- gamlss_mod_knots(pheno, knots=knots_list, sigma_knots=s_knots_list)
    }
    
    loop_count <- loop_count+1
    
    #if model isn't fit, skip to next loop
    if (is.null(model)) {
      message("model fitting failed, skipping to next iteration")
      next
    }

    saveRDS(model, file=paste0(save_path, "/model_objs/", pheno, "_mu", degree, "sig", sigma_degree, "_mod.rds"))
    
    #save centile fan plot
    print("creating centile fan plot")
    fan_plot <- centile_fan_resid(gamlssModel=model, df=df, x_var="logAge_days", color_var="sexMale",
                               get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                               sim_data_list = sim_df,
                               remove_cent_effect=c("fs_version", "study_site")) +
    ggtitle(paste(pheno, "\nsmoothed w/ mu.df=", degree, ", sigma.df=", sigma_degree))
    
    ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_mu", degree, "sig", sigma_degree, ".png"), fan_plot)
    
    #save worm plot
    print("creating worm plot")
    wp <- wp.taki(xvar=df$logAge_days, resid=resid(gamlssModel)) +
      ggtitle(paste(pheno, "\nsmoothed w/ mu.df=", degree, ", sigma.df=", sigma_degree))
    ggsave(file=paste0(save_path, "/worm_plots/", pheno, "_mu", degree, "sig", sigma_degree, ".png"), wp)
    
    #compile results
    sub_df <- cent_cdf(gamlssModel, df, "sexMale")
    sub_df$degrees <- paste0("mu_", degree, "sig_", sigma_degree)
    
    results_df <- rbind(results_df, sub_df)
  }
}

fwrite(results_df, file=paste0(save_path, "/cent_csvs/", pheno, "_results.csv"))

print(paste(sum(unique(results_df$degrees)), "of", loop_count, pheno, "models successful"))
