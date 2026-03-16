library(tidyr)
library(data.table)
library(dplyr)

#simplified version of centile_extreme_test.R to make sure denominators are calculated correctly

args <- commandArgs(trailingOnly = TRUE)
print(args)
dx_val <- as.character(args[1])
cv_sample <- as.character(args[2])

#base path
base_path <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
csv_path <- paste0(base_path, "cv_sample_", cv_sample, "_test/*totalFALSE*")

#get pheno list
lists <- list.files(paste0(base_path, "pheno_lists"), pattern = "\\.rds$", full.names = TRUE)
pheno_list <- do.call(c, lapply(lists, readRDS))

sum_df <- data.frame()
diffs_df <- data.frame()

#load patient and control centiles  CSVs
#custom fun to read and filt
fread_filt <- function(f){
  fread(f) %>%
    select(INDEX.ID, sex, dx_recode, contains("centile"))
}

#loop over pts and controls
for (grp in c("CN", "PT")){
  sum_df_grp <- data.frame()
  diffs_df_grp <- data.frame()
  denom_grp <- data.frame()
  
  #read in all centile csvs
  file_list <-list()
  
  for (pheno in pheno_list){
    f <- Sys.glob(paste0(csv_path, "/cent_csvs/", pheno, "_", grp, "_", dx_val, "_cent.csv"))
    if (length(f) == 0) {
      warning(paste("No file found for pheno:", pheno))
      next
    } else {
      file_list <- unique(unlist(c(file_list, f)))
    }}
  print(length(file_list))
  
  full_df <- rbindlist(lapply(file_list, fread_filt), fill=TRUE)
  
  cent_csv <- full_df %>%
    pivot_longer(
      cols = contains("centile"),
      names_to = c("pheno", "model"),
      names_pattern = "(.+)_centile_(full|null)"
    ) %>%
    filter(!is.na(value)) %>%
    mutate(ext = case_when(
      value < 0.05 ~ "low",
      value > 0.95 ~ "high",
      value >= 0.05 & value <= 0.95 ~ "mid",
      TRUE ~ NA))
  
  #summarise counts within each pheno
  sum_df_grp <- cent_csv %>%
    count(dx_recode, pheno, model, sex, ext)
  
  #summarise how many change across phenos
  diffs_df_grp <- cent_csv %>%
    group_by(INDEX.ID, pheno, sex, ext) %>%
    summarise(models = n_distinct(model), .groups = "drop") %>%
    filter(models == 1) %>% #subjects where just ONE model (sex-mod or sex-int) is extreme for that subj and pheno
    filter(ext != "mid") %>% #ignore when changes to/from mid so we're not double counting
    group_by(sex, ext) %>%
    summarise(n_change = n_distinct(INDEX.ID), .groups = "drop") #how many unique subjects have a change in at least one pheno?
  
  #get denominator - how many subjects per dx and sex were tested?  
  denom_grp <- full_df %>%
    select(INDEX.ID, sex, dx_recode) %>%
    distinct(INDEX.ID, .keep_all = TRUE) %>%
    count(dx_recode, sex) %>%
    rename(n_total = n)
  
  diffs_df_grp <- full_join(diffs_df_grp, denom_grp, by = c("sex")) %>%
    mutate(pct_change = (n_change / n_total) * 100)
  
  sum_df_grp <- full_join(sum_df_grp, denom_grp, by = c("dx_recode", "sex")) %>%
    mutate(pct = (n / n_total) * 100)
  
  sum_df <- rbind(sum_df, sum_df_grp)
  diffs_df <- rbind(diffs_df, diffs_df_grp)
}


#save
fwrite(sum_df, paste0(base_path, "cv_sample_", cv_sample, "_test/", dx_val, "_extcent_sum.csv"))
fwrite(diffs_df, paste0(base_path, "cv_sample_", cv_sample, "_test/", dx_val, "_extcent_diffs.csv"))