#Fit GAMLSS models to select from on CV samples

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)
library(perm)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))
options(warn = 1)

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to control csv
mod <- readRDS(args[2]) #test model fit on df
save_path <- as.character(args[3])

pt_df <- fread(paste0(base,"data/v3_pts_cleaned.csv"), stringsAsFactors = TRUE, na.strings = "") #path to patient data

mod$call$data <- "df"
mod$call$family <- mod$family[[1]]

##### PREP PATIENT DATAFRAME #####
pred_list <- list_predictors(mod)
pheno <- get_y(mod) %>% as.character()

all_list <- c(unlist(pred_list), pheno, "INDEX.ID", "dx", "dx_recode")
# print(length(all_list))
# print(all_list)

#drop extra variables
pt_df_clean <- pt_df %>%
  dplyr::select(all_of(all_list)) %>%
  na.omit()

##### CALCULATE CENTILES #####
print("calculating centiles...")
col_name <- paste0(pheno, "_cent")

#patients
pt_df_clean[[col_name]] <- pred_og_centile(mod, og.data=df, new.data=pt_df_clean)
print(range(pt_df_clean[[col_name]]))

#controls
df[[col_name]] <- pred_og_centile(mod, og.data=df)
print(range(df[[col_name]]))
#save
df$dx_recode <- "CN" #add dx_recode col to original df
df_full <- rbind(pt_df_clean, df, fill=TRUE)

fwrite(df_full, file=paste0(save_path, "/cent_csvs/", pheno, "_centiles.csv"))

##### LOOP OVER DX #####
print("testing disease effects...")
dx_list <- c("SCZ", "ALZ", "ASD", "MCI", "MDD", "GAD", "ADHD")

#make group col that accounts for dx and sex
# df_full <- df_full %>%
#   mutate(sex = ifelse(sexMale==0, "F", "M")) %>%
#   mutate(group = paste(sex, dx_recode, sep="_"))

# print(unique(df_full$group))

#holder
perm_df <- data.frame()

perm.opt <- permControl(nmc=1000,seed=123451, setSEED = TRUE)

for (dx in dx_list){
  print(paste("testing", dx))
  pt_cent <- pt_df_clean %>%
    mutate(dx_recode = as.character(dx_recode)) %>%
    filter(dx_recode == dx) %>%
    pull(col_name)
  
  #using code/opts from nature paper
  test <- permTS(df[[col_name]], pt_cent, paired=FALSE, alternative="two.sided", method ="exact.mc", control = perm.opt)

  test_df <- data.frame("dx" = dx,
                        "mean_diff" = test$estimate[[1]],
                        "p" = test$p.value,
                        "ci_low" = test$p.conf.int[[1]],
                        "ci_high" = test$p.conf.int[[2]],
                        "CN_median" = median(df[[col_name]]),
                        "pt_median" = median(pt_cent),
                        "pheno" = pheno
                        )
  perm_df <- rbind(perm_df, test_df)
}

fwrite(perm_df, file=paste0(save_path, "/cent_csvs/", pheno, "_cent_pt_test.csv"))

print("SUCCESS")
