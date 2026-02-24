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
df.og <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #df model fit on
mod_path <- as.character(args[2]) #test model fit on df
save_path <- as.character(args[3])
dx_val <- as.character(args[4])

#read in pt data
df_path <- paste0(base, "data/v3_pts_cleaned.csv")
df <- fread(df_path, stringsAsFactors = TRUE, na.strings = "") #df to fit centiles on


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

valid_sites <- unique(df.og$study_site)

#drop extra variables
df_clean <- df %>%
  dplyr::select(all_of(all_list)) %>%
  #remove extra sites
  dplyr::filter(study_site %in% valid_sites) %>%
  #filter to one dx
  filter(dx_recode == dx_val) %>%
  na.omit()

##### CALCULATE CENTILES #####
print("calculating centiles...")

df_cent <- lapply(names(mod_list), function(mn){
  m <- mod_list[[mn]]
  out_df <- pred_og_centile(
    m,
    og.data = df.og,
    new.data = df_clean,
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
df_full_cent <- bind_cols(df_clean, df_cent) %>%
  mutate(sex = ifelse(sexMale==0, "F", "M")) %>%
  #make sure base levels set at Female, Controls
  mutate(sex=factor(sex, 
                    levels = c("F", "M"),
                    ordered = TRUE))

#save centiles
fwrite(df_full_cent, file=paste0(save_path, "/cent_csvs/", pheno, "_all_", dx_val, "_cent.csv"))

##### DX #####
print("testing disease effects...")
dx_test_df <- data.frame()

#get sampling proportions for RESI
pi <- table(df_full_cent$sex)[["F"]]/nrow(df_full_cent)

for (mn in names(mod_list)){
  col_name <- paste(pheno, "std_score", mn, sep="_")
    
  #welch's t test - M vs F patients
  test_out <- t.test(df_full_cent[[col_name]] ~ df_full_cent$sex)
  test_d <- effsize::cohen.d(df_full_cent[[col_name]] ~ df_full_cent$sex)
  
  test_df <- data.frame("dx" = dx_val,
                        "model" = mn,
                        "pt.sex_tstat" = test_out$statistic,
                        "pt.sex_df" = test_out$parameter,
                        "pt.sex_p.val" = test_out$p.value,
                        "pt.sex_est" = test_out$estimate,
                        "pt.sex_ci_up" = test_out$conf.int[1],
                        "pt.sex_ci_low" = test_out$conf.int[1],
                        "pt.sex_d" = test_d$estimate,
                        "pt.sex_pi" = pi)
  
  dx_test_df <- rbind(dx_test_df, test_df)
  
}

fwrite(dx_test_df, file=paste0(save_path, "/cent_csvs/", pheno, "_cent_pt", dx_val, "_sex_test.csv"))

print("SUCCESS")
