#quick test for adding 1000 to phenos

set.seed(9999)

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlss.add)
library(mgcv)
library(gamlssTools)

#source("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/gamlss_fit_funs.R")

df_full <- fread("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/data/cv_sample_A.csv", stringsAsFactors=TRUE)
pheno_list <- readRDS("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/pheno_lists/global_vols.rds")


for (pheno in pheno_list) {
  pheno <- as.character(pheno)
  message("Processing: ", pheno)
  
  # Symbol for dplyr
  pheno_sym <- rlang::sym(pheno)
  
  # Select and mutate
  df_p <- df_full %>%
    dplyr::select(all_of(c(pheno, "fs_version_GM", "logAge_days", "sexMale_x_logAge", "sexMale", "study_site"))) %>%
    mutate(!!pheno_sym := (!!pheno_sym + 100000)) %>%
    na.omit()
  unscale <- function(x){x-100000}
  
  # Apply coverage truncation
  df_final <- trunc_coverage(df_p, vars = "logAge_days", n_min = 50)
  
  # Check df is valid
  stopifnot(is.data.frame(df_final))
  
  # Construct model formulas
  mu_formula <- as.formula(
    paste0(pheno, " ~ ga(~ s(sexMale_x_logAge, k=8, bs='ad') + ",
           "s(logAge_days, k=8, bs='ad') + sexMale + random(study_site), method = 'ML')")
  )
  
  sigma_formula <- ~ga(~s(sexMale_x_logAge, k=8, bs='ad') + 
                         s(logAge_days, k=8, bs='ad') + 
                         sexMale + random(study_site), method = 'ML')
  
  nu_formula <- ~ ga(~s(logAge_days, k=8, bs='ad') + 
                       sexMale, method = 'ML')
  
  print("fitting")
  
  # Fit the model
  test_mod <- gamlss(
    formula = mu_formula,
    sigma.formula = sigma_formula,
    nu.formula = nu_formula,
    family = "BCCG",
    data = df_final,
    control = gamlss.control(n.cyc = 400)
  )
  
  # Save model and centile plot
  saveRDS(test_mod, file = paste0("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/", pheno, "_test_modk8_nurefit2.rds"))
  
  fan_plot <- make_centile_fan(test_mod, 
                               df_final, 
                               x_var = "logAge_days", 
                               color_var = "sexMale", 
                               special_term = "sexMale_x_logAge = sexMale * logAge_days",
                               y_scale=unscale)
  
  ggsave(filename = paste0("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/", pheno, "_test_centilesk8_nurefit2.png"), fan_plot)
}
