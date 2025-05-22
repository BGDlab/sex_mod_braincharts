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
fs <- as.character(args[3])
total <- as.character(args[4])
save_path <- as.character(args[5])
log_pheno <- as.logical(args[6])
log_age <- as.logical(args[7])
sm <- as.character(args[8])

stopifnot(!is.null(total))

#loop over nu terms
if (log_age == TRUE & sm == "pb"){
  nu_list <- list(int = "1", 
                  site = "study_site", 
                  sex = "sexMale",
                  age = "logAge_days",
                  sexAge = "sexMale + logAge_days",
                  siteAge = "study_site + logAge_days", 
                  siteSex = "study_site + sexMale", 
                  siteAgeSex = "study_site + logAge_days + sexMale",
                  pbage = "pb(logAge_days, method='GAIC', k=log(nrow(df)), control = pb.control(order = 3))",
                  sex_pbAge = "sexMale + pb(logAge_days, method='GAIC', k=log(nrow(df)), control = pb.control(order = 3))",
                  site_pbAge = "study_site + pb(logAge_days, method='GAIC', k=log(nrow(df)), control = pb.control(order = 3))",
                  site_pbAgeSex = "study_site + pb(logAge_days, method='GAIC', k=log(nrow(df)), control = pb.control(order = 3)) + sexMale"
  )
  
} else if (log_age == TRUE & sm == "cs"){
  nu_list <- list(int = "1", 
                  site = "study_site", 
                  sex = "sexMale",
                  age = "logAge_days",
                  sexAge = "sexMale + logAge_days",
                  siteAge = "study_site + logAge_days", 
                  siteSex = "study_site + sexMale", 
                  siteAgeSex = "study_site + logAge_days + sexMale",
                  pbage = "cs(logAge_days)",
                  sex_pbAge = "sexMale + cs(logAge_days)",
                  site_pbAge = "study_site + cs(logAge_days)",
                  site_pbAgeSex = "study_site + cs(logAge_days) + sexMale"
  )
  
} else if (log_age == FALSE & sm == "pb"){
  nu_list <- list(int = "1", 
                  site = "study_site", 
                  sex = "sexMale",
                  age = "age_days",
                  sexAge = "sexMale + age_days",
                  siteAge = "study_site + age_days", 
                  siteSex = "study_site + sexMale", 
                  siteAgeSex = "study_site + age_days + sexMale",
                  pbage = "pb(age_days, method='GAIC', k=log(nrow(df)), control = pb.control(order = 3))",
                  sex_pbAge = "sexMale + pb(age_days, method='GAIC', k=log(nrow(df)), control = pb.control(order = 3))",
                  site_pbAge = "study_site + pb(age_days, method='GAIC', k=log(nrow(df)), control = pb.control(order = 3))",
                  site_pbAgeSex = "study_site + pb(age_days, method='GAIC', k=log(nrow(df)), control = pb.control(order = 3)) + sexMale"
  )
} else if (log_age == FALSE & sm == "cs"){
  nu_list <- list(int = "1", 
                  site = "study_site", 
                  sex = "sexMale",
                  age = "age_days",
                  sexAge = "sexMale + age_days",
                  siteAge = "study_site + age_days", 
                  siteSex = "study_site + sexMale", 
                  siteAgeSex = "study_site + age_days + sexMale",
                  pbage = "cs(age_days)",
                  sex_pbAge = "sexMale + cs(age_days)",
                  site_pbAge = "study_site + cs(age_days)",
                  site_pbAgeSex = "study_site + cs(age_days) + sexMale"
  )
}

#loop over fs moments
moment_list <- c("none", "mu", "both", "all")

#initialize empty lists
mod_list <- c()
results_df <- data.frame()
summary_df <- data.frame()

#FIT MODEL
for (fs_include in moment_list){
  
  for (total_include in moment_list){

    for (nu in nu_list){
    nu_name <- names(nu_list)[nu_list==nu]
    
    print(paste("fitting model with fs in", fs_include, 
                ", total in", total_include, 
                "and nu = ", nu_name))
    
    #FIT BASIC MODEL
    model <- gamlss_3lambda_etiv(pheno,
                            total_var=total, total_moment=total_include,
                            fs_ver=fs, fs_moment=fs_include, 
                            fam="BCCG",
                            nu_form=nu,
                            start.from = "mod_list[[1]]") #use first model as starting point
  
    #if model isn't fit, skip to next loop
    if (is.null(model)) {
      message("model fitting failed")
      next
    } else {
      message("model fit")
    }
    
    m_name <- paste(paste0("fs", fs_include), nu_name, paste0("total", total_include), sep="_")
    mod_list[[m_name]] <- model
    
    saveRDS(model, file=paste0(save_path, "/model_objs/", pheno, "_", m_name, "_mod.rds"))
  
   
    #COMPILE
      #BIC & AIC
      tmp_df <- data.frame(
        "AIC" = model$aic,
        "BIC" = model$sbc,
        "pheno" = pheno,
        "fs_moment" = fs_include,
        "total_moment" = total_include,
        "nu" = nu_name
      )
      summary_df <- rbind(summary_df, tmp_df)
    }
  }
}
expected <- length(nu_list)*length(moment_list)*length(moment_list)
print(paste(length(mod_list), "of", expected, "models fit"))

stopifnot(length(mod_list) > 0)

#SAVE CSVs
print("saving csvs")
fwrite(summary_df, file=paste0(save_path, "/model_sums/", pheno, "_summary.csv"))

print("finding lowest BIC")
best_bic <- summary_df %>%
  arrange(BIC) %>%
  slice_head(n=1)

best_bic$m_name <- paste(paste0("fs", best_bic$fs_moment), 
                         best_bic$nu, 
                         paste0("total", best_bic$total_moment), 
                         sep="_")

print(best_bic$m_name)

best_mod <- mod_list[[best_bic$m_name]]

#RENAME BEST MOD
file.rename(paste0(save_path, "/model_objs/", best_bic$pheno, "_", best_bic$m_name, "_mod.rds"),
            paste0(save_path, "/model_objs/", best_bic$pheno, "_", best_bic$m_name, "_BestMod.rds"))

print("compiling stats")

#CENTILE FAN PLOT
print("creating centile fan plots")
#sim data ONCE for centile fan plotting
print("simulate data for plotting")
sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")
sim_df2 <- sim_data(df, total, factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")

unscale <- function(x){10^x - 5}
fan_plot <- make_centile_fan(gamlssModel=best_mod, df=df, x_var="logAge_days", color_var="sexMale",
                             get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                             sim_data_list = sim_df,
                             remove_cent_effect="study_site",
                             remove_point_effect = c(total, "study_site"),
                             y_scale=unscale) +
  labs(title=paste(pheno, ",", best_bic$m_name),
     x ="log Age (days)",
     color = "Sex=Male", fill="Sex=Male")

ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_", best_bic$m_name, ".png"), fan_plot)

fan_plot <- make_centile_fan(gamlssModel=best_mod, df=df, x_var=total, color_var="sexMale",
                             get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                             sim_data_list = sim_df2,
			     remove_cent_effect="study_site",
                             remove_point_effect = c("logAge_days", "study_site"),
			     y_scale=unscale
                             ) +
  labs(title=paste(pheno, ",", best_bic$m_name),
       x = total,
       color = "Sex=Male", fill="Sex=Male")

ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_", best_bic$m_name, "total.png"), fan_plot)


#WORM PLOT
print("creating worm plot")
wp <- wp.taki(xvar=df$logAge_days, resid=resid(best_mod), n.inter=8)
ggsave(file=paste0(save_path, "/worm_plots/", pheno, "_", best_bic$m_name, ".png"), wp)

#centiles
results_df <- cent_cdf(best_mod, df, plot=FALSE, group="sexMale")
results_df$fs <- best_bic$fs_moment
results_df$total <- best_bic$total_moment
results_df$nu <- best_bic$nu

#centiles
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", pheno, "_", best_bic$m_name, "_results.csv"))

print("DONE")
