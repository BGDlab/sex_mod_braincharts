#Fit GAMLSS models to select from on CV samples

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))

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

family <- "BCTo"

stopifnot(total == "NULL")

#loop over nu terms
if (log_age == TRUE & sm == "pb"){
  nu_list <- list(
    int = "1",
    site = "study_site",
    sex = "sexMale",
    age = "logAge_days",
    sexAge = "sexMale + logAge_days",
    siteAge = "study_site + logAge_days",
    siteSex = "study_site + sexMale",
    siteAgeSex = "study_site + logAge_days + sexMale",
    pbAge = "pb(logAge_days, method='GAIC', k=log(nrow(df)))",
    sex_pbAge = "sexMale + pb(logAge_days, method='GAIC', k=log(nrow(df)))",
    site_pbAge = "study_site + pb(logAge_days, method='GAIC', k=log(nrow(df)))",
    site_pbAgeSex = "study_site + pb(logAge_days, method='GAIC', k=log(nrow(df))) + sexMale"
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
                  pbage = "pb(age_days, method='GAIC', k=log(nrow(df)))",
                  sex_pbAge = "sexMale + pb(age_days, method='GAIC', k=log(nrow(df)))",
                  site_pbAge = "study_site + pb(age_days, method='GAIC', k=log(nrow(df)))",
                  site_pbAgeSex = "study_site + pb(age_days, method='GAIC', k=log(nrow(df))) + sexMale"
  )
}


#loop over fs moments
fs_moment_list <- c("none", "mu", "both", "all")

#initialize empty lists
mod_count <- 0
first_mod <- NULL
results_df <- data.frame()
summary_df <- data.frame()

#FIT MODEL
for (fs_include in fs_moment_list){

  for (nu in nu_list){
  nu_name <- names(nu_list)[nu_list==nu]
  print(paste("fitting model with fs in", fs_include, "and nu = ", nu_name))
  
  m_name <- paste(fs_include, nu_name, family, sep="_")
  m_file <- paste0(save_path, "/model_objs/", pheno, "_", m_name, "_mod.rds")
  
  model <- NULL
  
  #CHECK IF MODEL EXISTS
  if (file.exists(m_file)){
    print("loading pre-fit model")
    model <- tryCatch({readRDS(m_file)
      }, error = function(e){
        message(e$message, "- trying again")
        tryCatch({readRDS(m_file)
          }, error = function(e){
            message(e$message, "- refit")
            NULL
            })
        })
    #double-check that loaded model did converge
    if (!is.null(model) && isFALSE(model$converged)) {
      print("loaded model was not converged")
      init_mod <- model
      print("fitting new gamlss model")
      #FIT BASIC MODEL
      if (sm == "pb" & log_age == TRUE){
        model <- gamlss_lambda(pheno,
                               fs_ver=fs,
                               fs_moment=fs_include,
                               fam=family,
                               nu_form=nu,
                               start.from = "init_mod") #use loaded model as starting point
        
      } else if (sm == "pb" & log_age == FALSE) {
        model <- gamlss_age(pheno,
                            fs_ver=fs,
                            fs_moment=fs_include,
                            fam=family,
                            nu_form=nu,
                            start.from = "init_mod") #use loaded model as starting point
        
      }
    }
  }
  
  #if no model, fit
  if (is.null(model)) {
    print("fitting new gamlss model")
    #FIT BASIC MODEL
    if (sm == "pb" & log_age == TRUE){
      model <- gamlss_lambda(pheno,
                             fs_ver=fs,
                             fs_moment=fs_include,
                             fam=family,
                             nu_form=nu,
                             start.from = "first_mod") #use first model as starting point
      
    } else if (sm == "pb" & log_age == FALSE) {
      model <- gamlss_age(pheno,
                          fs_ver=fs,
                          fs_moment=fs_include,
                          fam=family,
                          nu_form=nu,
                          start.from = "first_mod") #use first model as starting point
      
    } else if (sm == "cs" & log_age == TRUE){
      model <- gamlss_cs(pheno,
                         fs_ver=fs,
                         fs_moment=fs_include,
                         fam=family,
                         nu_form=nu,
                         start.from = "first_mod") #use first model as starting point
    } else if (sm == "cs" & log_age == FALSE){
      model <- gamlss_csage(pheno,
                            fs_ver=fs,
                            fs_moment=fs_include,
                            fam=family,
                            nu_form=nu,
                            start.from = "first_mod") #use first model as starting point
    }
  }
  
  #if model isn't fit, skip to next loop
  if (is.null(model)) {
    message("model fitting failed")
    next
  } else {
    message("model fit")
    mod_count <- mod_count + 1
  }
  saveRDS(model, file=m_file)
  
  #retain first successful model
  if (is.null(first_mod)){
    first_mod <- model
  }

  #COMPILE
    #BIC & AIC
    tmp_df <- data.frame(
      "AIC" = model$aic,
      "BIC" = model$sbc,
      "pheno" = pheno,
      "fs_moment" = fs_include,
      "nu" = nu_name
    )
    summary_df <- rbind(summary_df, tmp_df)
    
    #get back memory
    rm(model)
  }
}
expected <- length(nu_list)*length(fs_moment_list)
print(paste(mod_count, "of", expected, "models fit"))

stopifnot(mod_count > 0)

#SAVE CSVs
print("saving csvs")
fwrite(summary_df, file=paste0(save_path, "/model_sums/", pheno, "_summary.csv"))

print("finding lowest BIC")
best_bic <- summary_df %>%
  arrange(BIC) %>%
  slice_head(n=1) %>%
  tidyr::unite(m_name, c(fs_moment, nu)) %>%
  mutate(m_name = paste(m_name, family, sep="_"))

print(best_bic$m_name)

#RENAME BEST MOD
best_mod_file <- paste0(save_path, "/model_objs/", best_bic$pheno, "_", best_bic$m_name, "_BestMod.rds")
file.rename(paste0(save_path, "/model_objs/", best_bic$pheno, "_", best_bic$m_name, "_mod.rds"),
            best_mod_file)

print("compiling stats")

#READ BEST MOD BACK IN
best_mod <- readRDS(best_mod_file)
#re-write call info to be safe
best_mod$call$data <- "df"
best_mod$call$family <- family

#CENTILE FAN PLOT
#sim data ONCE for centile fan plotting
print("creating centile fan plot")
print("simulate data")
if (log_age == TRUE){
  sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")
  age_var <- "logAge_days"
} else {
  sim_df <- sim_data(df, "age_days", factor_var="sexMale", special_term = "sexMale_x_age = sexMale * age_days")
  age_var <- "age_days"
}

#if fs_version is included in BestMod, residualize from plot
if (fs %in% list_predictors(best_mod)){
  print("controlling for fs version")
  resid_terms <- c(fs, "study_site")
} else {
  resid_terms <- "study_site"
}

if (log_pheno==TRUE){
  unscale_fun <- unscale
} else {
  unscale_fun <- NULL
}

#plot
fan_plot <- make_centile_fan(gamlssModel=best_mod, df=df,
                             x_var=age_var,
                             color_var="sexMale",
                             get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                             sim_data_list = sim_df,
                             remove_point_effect = resid_terms,
                             y_scale=unscale_fun) +
  labs(title=paste(pheno, "\nsmoothed w/", sm, ",", best_bic$m_name),
       x = age_var,
       color = "Sex=Male", fill="Sex=Male")

ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_", best_bic$m_name, ".png"), fan_plot)

#WORM PLOT
print("creating worm plot")
wp <- wp.taki(xvar=df[[age_var]], resid=resid(best_mod), n.inter=6)$plot
wp <- wp + ggtitle(paste(pheno, "\nsmoothed w/", sm))
ggsave(file=paste0(save_path, "/worm_plots/", pheno, "_", best_bic$m_name, ".png"), wp)

#centiles
results_df <- cent_cdf(best_mod, df, plot=FALSE, group="sexMale")
results_df$fs <- best_bic$fs_moment
results_df$nu <- best_bic$nu

#centiles
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", pheno, "_", best_bic$m_name, "_results.csv"))

print("DONE")
