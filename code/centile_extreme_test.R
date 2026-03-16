library(tidyr)
library(data.table)
library(dplyr)

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

for (grp in c("CN", "PT")){
  all_diffs <- data.frame()
  denom_grp <- data.frame()
  
  for (pheno in pheno_list){
    f <- Sys.glob(paste0(csv_path, "/cent_csvs/", pheno, "_", grp, "_", dx_val, "_cent.csv"))
    if (length(f) == 0) {
      warning(paste("No file found for pheno:", pheno))
      next
    } else {
      print(paste(length(f), "pheno files found for", dx_val, grp))
    }
    
    cent_csv <- rbindlist(lapply(f, fread_filt)) %>%
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
    
    sum_df <- rbind(sum_df, cent_csv %>%
                      count(dx_recode, pheno, model, sex, ext))
    
    # collect subject-level disagreements across phenos
    all_diffs <- rbind(all_diffs, cent_csv %>%
                         group_by(INDEX.ID, dx_recode, sex, ext) %>%
                         summarise(models = n_distinct(model), .groups = "drop") %>%
                         filter(models == 1, ext != "mid"))
    
    denom_grp <- rbind(denom_grp, rbindlist(lapply(f, fread_filt)) %>%
                         distinct(INDEX.ID, dx_recode, sex))
  }
  
  # summarise changes once across all phenos
  csv_diffs_df <- all_diffs %>%
    group_by(dx_recode, sex, ext) %>%
    summarise(n_change = n_distinct(INDEX.ID), .groups = "drop")
  
  
  # finalize demographics
  denom_final <- denom_grp %>%
    distinct(INDEX.ID, dx_recode, sex) %>%
    count(dx_recode, sex) %>%
    rename(n_total=n)
  
  #add counts
  csv_diffs_df <- full_join(csv_diffs_df, denom_final, by = c("dx_recode", "sex")) %>%
    mutate(pct_change = (n_change / n_total) * 100)
  
  diffs_df <- rbind(diffs_df, csv_diffs_df)
  
  #add counts
  sum_df_final <- full_join(sum_df, denom_final, by = c("dx_recode", "sex"), 
                            relationship = "many-to-many") %>%
    mutate(pct = (n / n_total)*100)
  
}


#save
fwrite(sum_df_final, paste0(base_path, "cv_sample_", cv_sample, "_test/", dx_val, "_extcent_sum.csv"))
fwrite(diffs_df, paste0(base_path, "cv_sample_", cv_sample, "_test/", dx_val, "_extcent_diffs.csv"))
