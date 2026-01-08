#Fit GAMLSS models to select from on CV samples

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)
library(effsize)
library(broom)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))
options(warn = 1)

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to control csv
mod_path <- as.character(args[2]) #test model fit on df
save_path <- as.character(args[3])
dx_val <- as.character(args[4])

pt_df <- fread(paste0(base,"data/v3_pts_cleaned.csv"), stringsAsFactors = TRUE, na.strings = "") #path to patient data

##### READ IN MODS #####
mod <- readRDS(mod_path)
mod_list <- list(full = mod)

null_path <- sub("_full_", "_null_", mod_path)
if (file.exists(null_path)){
  null_mod <- readRDS(null_path)
  mod_list[["null"]] <- null_mod
}

null2_path <- sub("_full_", "_null2_", mod_path)
if (file.exists(null2_path)){
  null2_mod <- readRDS(null2_path)
  mod_list[["null2"]] <- null2_mod
}

## standardize model calls
mod_list <- lapply(mod_list, function(m) {
  m$call$data <- "df"
  m$call$family <- m$family[[1]]
  m
})

##### PREP PATIENT DATAFRAME #####
pred_list <- list_predictors(mod_list$full)
pheno <- get_y(mod_list$full) %>% as.character()

all_list <- c(unlist(pred_list), pheno, "INDEX.ID", "dx", "dx_recode")

# print(length(all_list))
# print(all_list)

#drop extra variables
pt_df_clean <- pt_df %>%
  dplyr::select(all_of(all_list)) %>%
  na.omit() %>%
  mutate(dx_recode = as.character(dx_recode)) %>%
  filter(dx_recode==dx_val)

print(paste(nrow(pt_df_clean), "patients"))
if(nrow(pt_df_clean)<1000){
  warning(paste(nrow(pt_df_clean), "patients, fewer than 1k"))
}
##### PREP CONTROL DATAFRAME #####
#find all sites with pt data
sites <- unique(pt_df_clean$study_site)
df_clean <- df %>%
  filter(study_site %in% sites)

print(paste(nrow(df_clean), "controls"))
if(nrow(df_clean)<1000){
  warning(paste(nrow(df_clean), "controls, fewer than 1k"))
}

##### JOIN #####
df_clean$dx_recode <- paste("CN", dx_val, sep="_") #add dx_recode col to original df
df_full <- rbind(pt_df_clean, df_clean, fill=TRUE)
print(unique(df_full$dx_recode))
stopifnot(length(unique(df_full$dx_recode))==2)

##### CALCULATE CENTILES #####
print("calculating centiles...")

df_cent <- lapply(names(mod_list), function(mn){
  m <- mod_list[[mn]]
  out_df <- pred_og_centile(
    m,
    og.data = df,
    new.data = df_full,
    get.std.scores = TRUE
  )
  # rename columns dynamically to include pheno and model type
  names(out_df) <- paste0(
    pheno, "_",
    names(out_df), "_",
    mn
  )
  out_df
})

#rejoin
df_full_cent <- bind_cols(df_full, df_cent) %>%
  mutate(sex = ifelse(sexMale==0, "F", "M"))

#save
fwrite(df_full_cent, file=paste0(save_path, "/cent_csvs/", pheno, "_", dx_val, "_cent.csv"))

##### DX #####
print("testing disease effects...")
dx_test_df <- data.frame()
dx_lm_df <- data.frame()

for (mn in names(mod_list)){
  col_name <- paste(pheno, "std_score", mn, sep="_")
    
  #welch's t test - patients vs controls
  test_out <- t.test(df_full_cent[[col_name]] ~ df_full_cent$dx_recode)
  test_d <- effsize::cohen.d(df_full_cent[[col_name]] ~ df_full_cent$dx_recode)
  
  #welch's t test - M vs F patients
  df_full_cent_pts <- df_full_cent %>%
    filter(dx_recode == dx_val)
  
  sex_test_out <- t.test(df_full_cent_pts[[col_name]] ~ df_full_cent_pts$sex)
  sex_test_d <- effsize::cohen.d(df_full_cent_pts[[col_name]] ~ df_full_cent_pts$sex)
  
  test_df <- data.frame("dx" = dx_val,
                        "model" = mn,
                        "case.control_tstat" = test_out$statistic,
                        "case.control_df" = test_out$parameter,
                        "case.control_p.val" = test_out$p.value,
                        "case.control_est" = test_out$estimate,
                        "case.control_ci_up" = test_out$conf.int[1],
                        "case.control_ci_low" = test_out$conf.int[1],
                        "case.control_d" = test_d$estimate,
                        
                        "pt.sex_tstat" = sex_test_out$statistic,
                        "pt.sex_df" = sex_test_out$parameter,
                        "pt.sex_p.val" = sex_test_out$p.value,
                        "pt.sex_est" = sex_test_out$estimate,
                        "pt.sex_ci_up" = sex_test_out$conf.int[1],
                        "pt.sex_ci_low" = sex_test_out$conf.int[1],
                        "pt.sex_d" = sex_test_d$estimate)
  
  dx_test_df <- rbind(dx_test_df, test_df)
  
  #lm
  lm_out <- lm(
    reformulate("dx_recode * sex", response = col_name),
    data = df_full_cent
  ) %>% tidy()
  
  lm_out$model <- mn
  
  dx_lm_df <- rbind(dx_lm_df, lm_out)
  
}

fwrite(dx_test_df, file=paste0(save_path, "/cent_csvs/", pheno, "_cent_pt", dx_val, "_test.csv"))
fwrite(dx_lm_df, file=paste0(save_path, "/cent_csvs/", pheno, "_cent_pt", dx_val, "_lm.csv"))

print("SUCCESS")
