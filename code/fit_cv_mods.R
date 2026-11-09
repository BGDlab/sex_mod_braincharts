#Fit GAMLSS models to select from on CV samples

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)
library(parallel)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/"
source(paste0(base, "code/gamlss_fit_funs.R"))

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- droplevels(fread(args[1], stringsAsFactors = TRUE, na.strings = c("", '""'))) #path to csv; drop unused factor levels
pheno <- as.character(args[2])
fs <- as.character(args[3])
total <- as.character(args[4])
save_path <- as.character(args[5])
log_age <- as.logical(args[6])
sm <- as.character(args[7])
#cores to fit the model grid on (pass $SLURM_CPUS_PER_TASK; singularity --cleanenv hides it from R); default 1 = serial
n_cores <- if (length(args) >= 8 && !is.na(suppressWarnings(as.integer(args[8])))) as.integer(args[8]) else 1L

family <- "BCTo"

stopifnot(total == "NULL")

#SEX CODING: pb(logAge_days) + pvc(logAge_days, by = sexEffect) (see gamlss_fit_funs.R); only set up for log-age pb models
stopifnot(log_age == TRUE, sm == "pb")
stopifnot("df has no sexEffect column (rerun lbcc_eda.Rmd)" = SEX_VAR %in% names(df))

#loop over nu terms
if (log_age == TRUE & sm == "pb"){
  nu_list <- list(
    int = "1",
    site = "study_site",
    sex = "sexEffect",
    age = "logAge_days",
    sexAge = "sexEffect + logAge_days",
    siteAge = "study_site + logAge_days",
    siteSex = "study_site + sexEffect",
    siteAgeSex = "study_site + logAge_days + sexEffect",
    pbAge = "pb(logAge_days, method='GAIC', k=log(nrow(df)))",
    sex_pbAge = "sexEffect + pb(logAge_days, method='GAIC', k=log(nrow(df)))",
    site_pbAge = "study_site + pb(logAge_days, method='GAIC', k=log(nrow(df)))",
    site_pbAgeSex = "study_site + pb(logAge_days, method='GAIC', k=log(nrow(df))) + sexEffect"
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

#model grid, in the order the old nested loops fit it: fs moment, then nu (fastest)
grid <- expand.grid(nu_name = names(nu_list), fs_include = fs_moment_list,
                    stringsAsFactors = FALSE)[, c("fs_include", "nu_name")]

#fit (or load) grid model i, save it, and return its summary row (NULL summary if it failed to fit)
#return_model = TRUE also returns the model (only used for first_mod: workers don't ship models back)
fit_grid_mod <- function(i, return_model = FALSE){
  fs_include <- grid$fs_include[i]
  nu_name <- grid$nu_name[i]
  nu <- nu_list[[nu_name]]
  print(paste("fitting model with fs in", fs_include, "and nu = ", nu_name))
  
  m_name <- paste(fs_include, nu_name, family, sep="_")
  m_file <- paste0(save_path, "/model_objs/", pheno, "_", m_name, "_mod.rds")
  
  model <- NULL
  
  #CHECK IF MODEL EXISTS
  if (file.exists(m_file)){
    print(paste(m_name, ": loading pre-fit model"))
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
      print(paste(m_name, ": loaded model was not converged"))
      #gamlss_lambda() evaluates its call text with lookup through the global env, so init_mod has to be
      #global (as in the old top-level loop); in a forked worker this only touches that worker's copy
      assign("init_mod", model, envir = globalenv())
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
      
    }
  }
  
  #if model isn't fit, skip to next model
  if (is.null(model)) {
    message(m_name, ": model fitting failed")
    return(list(summary = NULL, model = NULL))
  }
  message(m_name, ": model fit")
  saveRDS(model, file=m_file)
  
  #COMPILE
  #BIC & AIC
  tmp_df <- data.frame(
    "AIC" = model$aic,
    "BIC" = model$sbc,
    "pheno" = pheno,
    "fs_moment" = fs_include,
    "nu" = nu_name,
    "mu_pvc_lambda" = pvc_lambda(model, "mu"), #1e7 = capped, i.e. linear sex difference in age
    "sigma_pvc_lambda" = pvc_lambda(model, "sigma")
  ) %>%
    cbind(df_demo_stats(df)) #age (days) & sex of fit data
  
  list(summary = tmp_df, model = if (return_model) model)
}

#FIT MODELS
fits <- vector("list", nrow(grid))
first_mod <- NULL

#fit in grid order until one model fits: that first successful model is the starting point for every later fit
i <- 0
while (is.null(first_mod) && i < nrow(grid)){
  i <- i + 1
  fits[[i]] <- fit_grid_mod(i, return_model = TRUE)
  first_mod <- fits[[i]]$model
  fits[[i]]$model <- NULL
}

#the remaining fits only depend on first_mod, so run them in parallel (forked workers inherit df, first_mod
#and the patched gamlss namespace); mc.preschedule = FALSE hands out one model at a time, since fit times vary a lot
rest <- seq_len(nrow(grid))[-seq_len(i)]
if (length(rest) > 0){
  print(paste("fitting", length(rest), "remaining models on", n_cores, "cores"))
  fits[rest] <- mclapply(rest, fit_grid_mod, mc.cores = n_cores, mc.preschedule = FALSE)
}

#a worker that hit an uncaught error (which stopped the old serial loop) returns a try-error, and one that was
#killed (e.g. out of memory) returns NULL: stop before picking a best model from an incomplete set
#(models that did fit are saved, so a rerun only refits the missing ones)
bad <- which(vapply(fits, function(f) is.null(f) || inherits(f, "try-error"), logical(1)))
if (length(bad) > 0){
  for (b in bad){
    message("model ", paste(grid[b, ], collapse = "/"), " errored: ",
            if (is.null(fits[[b]])) "worker died (out of memory?)" else as.character(fits[[b]]))
  }
  stop(length(bad), " model fit(s) errored; rerun to refit them")
}

summary_df <- do.call(rbind, lapply(fits, `[[`, "summary"))
mod_count <- if (is.null(summary_df)) 0 else nrow(summary_df)
expected <- nrow(grid)
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
  sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = SEX_SPECIAL_TERM)
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

#plot
fan_plot <- make_centile_fan(gamlssModel=best_mod, df=df,
                             x_var=age_var,
                             color_var="sexMale",
                             get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                             sim_data_list = sim_df,
                             remove_point_effect = resid_terms) +
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
