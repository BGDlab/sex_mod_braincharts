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
base_mod <- readRDS(args[2])
save_path <- as.character(args[3])
total <- as.character(args[4])

filename_no_ext <- sub("\\.[^.]*$", "", basename(args[2]))
filename <- sub("train", "test", filename_no_ext)
file_full <- paste0(save_path, "/model_objs/", filename, "_full_mod.rds")

# #check if this pheno is already run, and if so, end
# if (file.exists(paste0(save_path, "/model_sums/", filename, "_LRtest.csv"))){
#   stop("Already tested, skipping pheno")
# }

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
    dplyr::select(all_of(c(pred_list, pheno, "weight"))) %>%
    na.omit() %>%
    trunc_coverage(age_var) #drop points at ends if too sparse
} else {
  df <- full_df %>%
    dplyr::select(all_of(c(pred_list, pheno, "weight"))) %>%
  na.omit() %>%
    trunc_coverage(c(total, age_var)) #drop points at ends if too sparse
}

##### FIT TEST MODEL #####
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
  model <- gamlss_lambda_rep(base_mod, 
                             null_mod="false",
                             keep_lambdas=TRUE,
                             weight="weight")
}

#if model isn't fit, skip to next loop
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

model$call$data <- df
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
    labs(title=paste(pheno, "QC-weighted validation model"),
         x =age_var,
         color = "Sex=Male", fill="Sex=Male")
    
    ggsave(file=paste0(save_path, "/centile_plots/", filename, ".png"), fan_plot)
    
  #WORM PLOT
    print("creating worm plot")
    wp <- wp.taki(xvar=df[[age_var]], resid=resid(model), n.inter=6)$plot +
      ggtitle(paste(pheno, "QC-weighted validation model"))
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
    )

#SAVE CSVs
print("saving csvs")

#centiles
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", filename, "_centiles.csv"))

#BIC & AIC
fwrite(summary_df, file=paste0(save_path, "/model_sums/", filename, "_summary.csv"))

##################
#FIT NULL MODEL
print("fitting null model")
null_model <- gamlss_lambda_rep(base_mod, null_mod="true", weight="weight")
#test saving null 
file_null <- paste0(save_path, "/model_objs/", filename, "_null_mod.rds")
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
if (total != FALSE) {
  print("fitting null model of all sex effects")
  null_model2 <- gamlss_lambda_rep(base_mod, null_mod="allSex", weight="weight")
  
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

print("SUCCESS")
