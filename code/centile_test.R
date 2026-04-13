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
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #df to fit centiles on
df.og <- fread(args[2], stringsAsFactors = TRUE, na.strings = "") #df model fit on
mod_path <- as.character(args[3]) #test model fit on df
save_path <- as.character(args[4])
dx_val <- as.character(args[5])

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
df_clean <- df %>%
  dplyr::select(all_of(all_list)) %>%
  na.omit()
stopifnot(length(unique(df_clean$dx_recode))==2)

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

cn_level <- grep("^CN_", levels(df_full_cent$dx_recode), value = TRUE)
pt_level <- setdiff(levels(df_full_cent$dx_recode), cn_level)

# Reorder with CN_ as base
df_full_cent <- df_full_cent %>%
  mutate(dx_recode = factor(dx_recode, 
                            levels = c(cn_level, pt_level),
                            ordered = TRUE))

#save patient and control centiles separately, otherwise files are too big to read in
df_cn <- df_full_cent %>%
  filter(dx_recode==cn_level)
fwrite(df_cn, file=paste0(save_path, "/cent_csvs/", pheno, "_CN_", dx_val, "_cent.csv"))

df_pt <- df_full_cent %>%
  filter(dx_recode!=cn_level)
fwrite(df_pt, file=paste0(save_path, "/cent_csvs/", pheno, "_PT_", dx_val, "_cent.csv"))