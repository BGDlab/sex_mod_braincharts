#Fit GAMLSS models to select from on CV samples

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)
library(ggplot2)
library(ggnewscale)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/"
source(paste0(base, "code/gamlss_fit_funs.R"))
source(paste0(base, "code/sex_diff_funs.R"))
options(warn = 1)

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- droplevels(fread(args[1], stringsAsFactors = TRUE, na.strings = c("", '""'))) #path to csv; drop unused factor levels
base_mod <- readRDS(args[2])
save_path <- as.character(args[3])
total <- as.logical(args[4])
sim_df_uniform <- readRDS(args[5])

filename_no_ext <- sub("\\.[^.]*$", "", basename(args[2]))
filename <- sub("BestMod", "test", filename_no_ext)
file_full <- paste0(save_path, "/model_objs/", filename, "_full_mod.rds")

base_mod$call$data <- "df"
base_mod$call$family <- base_mod$family[[1]]

#FIT BASIC MODEL
#check if this pheno is already run, and if so, load
model <- NULL

if (file.exists(file_full)){
  print("loading pre-fit model")
  model <- tryCatch({readRDS(file_full)
  }, error = function(e){
    message(e$message, "- trying again")
    tryCatch({readRDS(file_full)
    }, error = function(e){
      message(e$message, "- refit")
      NULL
    })
  })
}
if (is.null(model)){
  model <- gamlss_lambda_rep(base_mod, null_mod="false", n.cyc=1800)
}

#if model doesn't fit, fail
if (is.null(model)) {
  message("model fitting failed")
  stop("model fitting failed")
}

  #compare lambdas
  print("original lambdas:")
  print(base_mod$mu.lambda)
  print(base_mod$sigma.lambda)
  print(base_mod$nu.lambda)  

  print("refit lambdas:")
  print(model$mu.lambda)
  print(model$sigma.lambda)
  print(model$nu.lambda)
  
  mu_diff <- base_mod$mu.lambda - model$mu.lambda
  sig_diff <- base_mod$sigma.lambda - model$sigma.lambda
  stopifnot(c(mu_diff, sig_diff) < 0.01)
  if (!is.null(base_mod$nu.lambda)){
	nu_diff <- base_mod$nu.lambda - model$nu.lambda
	stopifnot(nu_diff < 0.01)}

print (paste("saving to", file_full))
saveRDS(model, file=file_full)

model$call$data <- "df"
model$call$family <- model$family[[1]]

#CENTILE FAN PLOT
pheno <- model$mu.terms[[2]] %>% as.character()

pred_list <- list_predictors(model)

#check if age is log-scaled
if ("logAge_days" %in% pred_list){
  age_var <- "logAge_days"
  print("simulate data for plotting")
  sim_df <- sim_data(df, age_var, factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")
  vars_of_interest <- c("sexMale_x_logAge", age_var, "sexMale")
} else {
  print("simulate data for plotting")
  age_var <- "age_days"
  sim_df <- sim_data(df, age_var, factor_var="sexMale", special_term = "sexMale_x_age = sexMale * age_days")
  vars_of_interest <- c("sexMale_x_age", age_var, "sexMale")
}

#define nuisance covars to residualize from points
resid_terms <- setdiff(pred_list, vars_of_interest)

#get split-half (A/B) from save_path, e.g. .../cv_sample_A_test/...
split <- sub(".*cv_sample_([AB])_test.*", "\\1", save_path)
stopifnot(split %in% c("A", "B"))

birth <- log(280, base=10)
pt_show <- !isTRUE(total)
title_text <- paste0(pheno," test model in Split-Half ", split)
print("creating centile fan plot")
  fan_plot <- make_centile_fan(gamlssModel=model, 
                               df=df, 
                               x_var=age_var, 
                               color_var="sexMale",
                               get_peaks=TRUE, 
                               desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                               show_points=pt_show,
                               sim_data_list = sim_df,
                               label_centiles = "legend",
                               point_color_manual = c('0' = "#F4A15BFF", '1' = "#8CB3D1FF"),
                               color_manual = c("0" = "#c05f0d", "1" = "#38688D"),
                               remove_point_effect = resid_terms,
                               x_axis="log_lifespan_fetal") +
    theme_linedraw() +
    scale_color_manual(name = "Sex",
                       values = c('0' = "#c05f0d", '1' = "#38688D"),
                       labels = c("Female", "Male")) +
    guides(fill = "none",
           fill_ggnewscale_1 = "none", colour_ggnewscale_1 = "none") +
    geom_vline(xintercept=birth, color="gray20", linetype = "dashed") +
    labs(title= title_text,
         x = "Age (years)")
    
    ggsave(file=paste0(save_path, "/centile_plots/", filename, ".png"), fan_plot, width = 8, height = 5)
    
  #WORM PLOT
    print("creating worm plot")
    wp <- wp.taki(xvar=df[[age_var]], resid=resid(model), n.inter=6)$plot +
      ggtitle(paste(pheno, "validation model"))
    ggsave(file=paste0(save_path, "/worm_plots/", filename, ".png"), wp)
    
  #COMPILE
    print("compiling stats")
    #centiles
    results_df <- cent_cdf(model, df, plot=FALSE, group="sexMale")
    
    #BIC & AIC
    summary_df <- data.frame(
      "AIC" = model$aic,
      "BIC" = model$sbc,
      "pheno" = pheno
    ) %>%
      cbind(df_demo_stats(df)) #age (days) & sex of fit data

#SAVE CSVs
print("saving csvs")

#centiles
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", filename, "_centiles.csv"))

#BIC & AIC
fwrite(summary_df, file=paste0(save_path, "/model_sums/", filename, "_summary.csv"))

##################
#FIT NULL MODEL
print("fitting null model")

#check if null is already run, and if so, load
null_model <- NULL
file_null <- paste0(save_path, "/model_objs/", filename, "_null_mod.rds")

if (file.exists(file_null)){
  print("loading pre-fit model")
  null_model <- tryCatch({readRDS(file_null)
  }, error = function(e){
    message(e$message, "- trying again")
    tryCatch({readRDS(file_null)
    }, error = function(e){
      message(e$message, "- refit")
      NULL
    })
  })
}
if (is.null(null_model)){
  null_model <- gamlss_lambda_rep(base_mod, null_mod="true", n.cyc=1200)
}

#if model doesn't fit, fail
if (is.null(null_model)) {
  message("null model fitting failed")
  stop("null model fitting failed")
}

print (paste("saving to", file_null))
saveRDS(null_model, file=file_null)

test_out <- LR.test(null_model, model, print=FALSE) #significance test
f2 <- cohens_f2_local(model, null_model) #effect size

#TEST
test_df <- data.frame(
  "chi" = test_out$chi,
  "df" = test_out$df,
  "p_val" = test_out$p.val,
  "fsq" = f2,
  "pheno" = pheno,
  "effect" = "sex_age"
)

##################
#FIT NULL MODEL WITH NO SEX-EFFECT FOR TBV-CORRECTED MODELS
if (total == TRUE) {
  print("fitting null model of all sex effects")
  null_model2 <- gamlss_lambda_rep(base_mod, null_mod="allSex", n.cyc=1200)
  
  file_null2 <- paste0(save_path, "/model_objs/", filename, "_null2_mod.rds")
  print (paste("saving to", file_null2))
  saveRDS(null_model2, file=file_null2)
  
  test_out2 <- LR.test(null_model2, model, print=FALSE) #significance test
  f2_2 <- cohens_f2_local(model, null_model2) #effect size
  
  #TEST
  test_df2 <- data.frame(
    "chi" = test_out2$chi,
    "df" = test_out2$df,
    "p_val" = test_out2$p.val,
    "fsq" = f2_2,
    "pheno" = pheno,
    "effect" = "sex_all"
  )
  
  test_df <- rbind(test_df, test_df2)
}

fwrite(test_df, file=paste0(save_path, "/model_sums/", filename, "_LRtest.csv"))

####################
#GET SEX-DIFFERENCES
# Predict 50th cent & sigma values for each simulated level of factor_var


# Process sim_df_uniform (shared grid for cross-pheno/split comparison)
result_df <- process_sex_diffs(sim_data_list = sim_df_uniform, model = model, df = df,
                               age_var = age_var, pred_list = pred_list, pheno = pheno)
# Save results
print("saving uniform sex diffs")
fwrite(result_df, file=paste0(save_path, "/cent_csvs/", pheno, "_uniform_sexdiffs.csv"))

print("SUCCESS")
