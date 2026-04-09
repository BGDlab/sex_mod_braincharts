library(tidyr)
library(data.table)
library(dplyr)
library(purrr)
library(EnvStats)

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
  
  # --- CHUNKED pivot_longer to avoid vec_interleave_indices() size limit ---
  chunk_size <- 50  # number of files per chunk; tune down if still hitting the limit
  
  cent_cols <- grep("centile", names(full_df), value = TRUE)
  id_cols   <- c("INDEX.ID", "sex", "dx_recode")
  phenos_in_data <- unique(sub("_centile_(full|null)$", "", cent_cols))
  
  chunks <- split(phenos_in_data, ceiling(seq_along(phenos_in_data) / chunk_size))
  
  cent_csv <- rbindlist(lapply(chunks, function(pheno_chunk) {
    cols_chunk <- cent_cols[sub("_centile_(full|null)$", "", cent_cols) %in% pheno_chunk]
    sub_df <- full_df[, c(id_cols, cols_chunk), with = FALSE]
    
    sub_df %>%
      pivot_longer(
        cols = all_of(cols_chunk),
        names_to = c("pheno", "model"),
        names_pattern = "(.+)_centile_(full|null)"
      ) %>%
      filter(!is.na(value)) %>%
      mutate(ext = case_when(
        value < 0.05  ~ "low",
        value > 0.95  ~ "high",
        value >= 0.05 & value <= 0.95 ~ "mid",
        TRUE ~ NA
      ))
  }), fill = TRUE)
  # -------------------------------------------------------------------------
  
  #summarize counts within each pheno
  sum_df_grp <- cent_csv %>%
    count(dx_recode, pheno, model, sex, ext)
  
  #summarize
  diffs_df_grp <- cent_csv %>%
    group_by(INDEX.ID, pheno, sex, ext) %>%
    summarise(models = n_distinct(model), .groups = "drop") %>%
    filter(models == 1) %>% #subjects where just ONE model (sex-mod or sex-int) is extreme for that subj and pheno
    filter(ext != "mid") %>% #ignore when changes to/from mid so we're not double counting
    #how many phenos changed for a given subject?
    group_by(INDEX.ID, sex) %>%
    summarise(n_pheno_change = n_distinct(pheno), .groups = "drop") %>%
    group_by(sex) %>%
    summarise(n_change = n_distinct(INDEX.ID), #how many unique subjects have a change in at least one pheno?
              n_pheno_change_total = sum(n_pheno_change), #how many phenos change at all?
              #of subjects who changed in at least one pheno, how many phenos changed?
              min_phenos = min(n_pheno_change),
              max_phenos = max(n_pheno_change),
              mean_phenos = mean(n_pheno_change),
              .groups = "drop") 
  
  #get denominator - how many subjects per dx and sex were tested?  
  denom_grp <- full_df %>%
    select(INDEX.ID, sex, dx_recode) %>%
    distinct(INDEX.ID, .keep_all = TRUE) %>%
    count(dx_recode, sex) %>%
    rename(n_total = n)
  
  diffs_df_grp <- full_join(diffs_df_grp, denom_grp, by = c("sex")) %>%
    mutate(pct_change = (n_change / n_total) * 100,
           pct_pheno_change = (n_pheno_change_total/(n_total*length(pheno_list)))*100)
  
  sum_df_grp <- full_join(sum_df_grp, denom_grp, by = c("dx_recode", "sex")) %>%
    mutate(pct = (n / n_total) * 100,
           dx = grp)
  
  sum_df <- rbind(sum_df, sum_df_grp)
  diffs_df <- rbind(diffs_df, diffs_df_grp)
}


#save
sum_df$dx_tested <- dx_val
sum_df$cv_sample <- cv_sample
diffs_df$dx_tested <- dx_val
diffs_df$cv_sample <- cv_sample

fwrite(sum_df, paste0(base_path, "cv_sample_", cv_sample, "_test/", dx_val, "_extcent_sum.csv"))
fwrite(diffs_df, paste0(base_path, "cv_sample_", cv_sample, "_test/", dx_val, "_extcent_diffs.csv"))

#SIGNIFICANCE TESTING using Fisher's Exact test of proportions
run_fisher <- function(n_CN, n_PT, n_total_CN, n_total_PT) {
  twoSamplePermutationTestProportion(
    x = c(n_CN, n_PT),
    y = c(n_total_CN, n_total_PT),
    x.and.y = "Number Successes and Trials",
    alternative = "two.sided"
  )
}

#prep df
fisher_df <- sum_df %>%
  filter(ext != "mid") %>%
  pivot_wider(id_cols=c(pheno, model, sex, ext), names_from=dx, values_from=c(n, n_total), values_fill = 0) %>%
  #collapsing across sex for now
  group_by(pheno, model, ext) %>%
  summarise(n_CN = sum(n_CN),
            n_PT = sum(n_PT),
            n_total_CN = sum(n_total_CN),
            n_total_PT = sum(n_total_PT)) 
#run stats
fisher_results <- pmap(
  list(fisher_df$n_CN, fisher_df$n_PT, fisher_df$n_total_CN, fisher_df$n_total_PT),
  run_fisher
)

fisher_final <- fisher_df %>%
  ungroup() %>%
  mutate(
    fisher_stat = map_dbl(fisher_results, ~ .x$statistic),
    fisher_p    = map_dbl(fisher_results, ~ .x$p.value)
  )

fisher_final$dx_tested <- dx_val
fisher_final$cv_sample <- cv_sample

fwrite(fisher_final, paste0(base_path, "cv_sample_", cv_sample, "_test/", dx_val, "_extcent_fisher.csv"))
