library(tidyr)
library(data.table)
library(dplyr)
library(purrr)
library(EnvStats)

#simplified version of centile_extreme_test.R to make sure denominators are calculated correctly

args <- commandArgs(trailingOnly = TRUE)
print(args)
dx_val <- as.character(args[1])
#cv_sample <- as.character(args[2])

#set paths
base_path <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
csv_path <- paste0(base_path, "cv_sample_?_test/*totalFALSE*")
save_path <- paste0(base_path, "dx_tests/")

#get pheno list
lists <- list.files(paste0(base_path, "pheno_lists"), pattern = "\\.rds$", full.names = TRUE)
pheno_list <- do.call(c, lapply(lists, readRDS))

#READ IN AND AVERAGE PT CENTILES
fread_filt <- function(f, string){
  fread(f) %>%
    select(INDEX.ID, sexMale, dx_recode, logAge_days, matches(string)) %>%
    mutate(sex = ifelse(sexMale==0, "F", "M"))
}

pt_df_list <- c()
#average within pheno
for (pheno in pheno_list) {
  f_list <- Sys.glob(paste0(csv_path, "/cent_csvs/", pheno, "_PT_", dx_val, "_cent.csv"))
  if (length(f_list) != 2) {
    warning(paste("2 files not found for pheno:", pheno))
    next
  } else {
    pheno_df <-  rbindlist(lapply(f_list, fread_filt, "std_score"), fill=TRUE)
  }
  
  pheno_mean <- pheno_df %>%
    group_by(INDEX.ID, sex, dx_recode) %>%
    summarise(across(contains("std_score"), mean), .groups = "drop") %>%
    # convert mean std scores back to centiles
    mutate(
      across(
        .cols = contains("std_score") & 
          !ends_with("_diff") & 
          !ends_with("_diff2"), #don't back-convert diff cols
        .fns = pnorm,
        .names = "{gsub('std_score', 'centile', .col)}"
      )
    ) %>%
    ungroup()
  
  pt_df_list <- c(pt_df_list, list(pheno_mean))
}
print(length(pt_df_list))

pt_df <- pt_df_list %>% purrr::reduce(dplyr::full_join, by = c("INDEX.ID", "sex", "dx_recode"))
stopifnot(nrow(pt_df) == length(unique(pt_df$INDEX.ID)))
fwrite(pt_df, paste0(save_path, dx_val, "_mean_scores.csv")) #figure out where to save

cn_df_list <- c()
##READ IN CONTROLS
for (pheno in pheno_list) {
  f_list <- Sys.glob(paste0(csv_path, "/cent_csvs/", pheno, "_CN_", dx_val, "_cent.csv"))
  if (length(f_list) != 2) {
    warning(paste("2 files not found for pheno:", pheno))
    next
  } else {
    pheno_df <-  rbindlist(lapply(f_list, fread_filt, "std_score|centile"), fill=TRUE)
  }
  
  cn_df_list <- c(cn_df_list, list(pheno_df))
}
print(length(cn_df_list))

cn_df <- cn_df_list %>% purrr::reduce(dplyr::full_join, by = c("INDEX.ID", "sex", "dx_recode"))
stopifnot(nrow(cn_df) == length(unique(cn_df$INDEX.ID)))

#T TEST TO COMPARE MEAN Z SCORES IN EACH PHENO
#also tests differences in between-model changes
std_cols <- names(cn_df)[grepl("std_score", names(cn_df))]

t_results <- map_dfr(std_cols, function(col) {
  #print(col)
  x <- na.omit(cn_df[[col]])
  y <- na.omit(pt_df[[col]])
  
  if (length(x) == 0 || length(y) == 0) {
  message("Skipping ", col)
  return(NULL)
  }  

  test_out <- t.test(x, y)
  d_out <- effsize::cohen.d(x, y)
  
  tibble(
    variable = col,
    t_stat = test_out$statistic,
    p_value = test_out$p.value,
    mean_df1 = mean(x, na.rm = TRUE),
    mean_df2 = mean(y, na.rm = TRUE),
    n_cn = sum(!is.na(x)),
    n_pt = sum(!is.na(y)),
    d = d_out$estimate
  )
})

#save
fwrite(t_results, file=paste0(save_path, dx_val, "_casecontrol_test.csv"))

#TEST EXTREMENESS & CHANGES IN EXTREMENESS ACROSS MODELS IN CASES V CONTROLS
sum_df <- data.frame()
diffs_df <- data.frame()
#loop over pts and controls
for (df in list(cn_df, pt_df)){
  sum_df_grp <- data.frame()
  diffs_df_grp <- data.frame()
  denom_grp <- data.frame()
  
  #get dx
  x <- unique(df$dx_recode)
  dx_grp <- ifelse(grepl("CN", x), "CN", "PT")
  
  # --- CHUNKED pivot_longer to avoid vec_interleave_indices() size limit ---
  chunk_size <- 50  # number of files per chunk; tune down if still hitting the limit
  
  cent_cols <- grep("centile", names(df), value = TRUE)
  id_cols   <- c("INDEX.ID", "sex", "dx_recode")
  phenos_in_data <- unique(sub("_centile_(full|null)$", "", cent_cols))
  
  chunks <- split(phenos_in_data, ceiling(seq_along(phenos_in_data) / chunk_size))
  
  cent_csv <- rbindlist(lapply(chunks, function(pheno_chunk) {
    cols_chunk <- cent_cols[sub("_centile_(full|null)$", "", cent_cols) %in% pheno_chunk]
    sub_df <- df[, c(id_cols, cols_chunk), with = FALSE]
    
    sub_df %>%
      pivot_longer(
        cols = all_of(cols_chunk),
        names_to = c("pheno", "model"),
        names_pattern = "(.+)_centile_(full|null)"
      ) %>%
      filter(!is.na(value) & !is.na(model)) %>%
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
  denom_grp <- df %>%
    select(INDEX.ID, sex, dx_recode) %>%
    distinct(INDEX.ID, .keep_all = TRUE) %>%
    count(dx_recode, sex) %>%
    rename(n_total = n)
  
  diffs_df_grp <- full_join(diffs_df_grp, denom_grp, by = c("sex")) %>%
    mutate(pct_change = (n_change / n_total) * 100,
           pct_pheno_change = (n_pheno_change_total/(n_total*length(pheno_list)))*100)
  
  sum_df_grp <- full_join(sum_df_grp, denom_grp, by = c("dx_recode", "sex")) %>%
    mutate(pct = (n / n_total) * 100,
           dx = dx_grp)
  
  sum_df <- rbind(sum_df, sum_df_grp)
  diffs_df <- rbind(diffs_df, diffs_df_grp)
}


#save
sum_df$dx_tested <- dx_val
diffs_df$dx_tested <- dx_val

fwrite(sum_df, file=paste0(save_path, dx_val, "_extcent_sum.csv"))
fwrite(diffs_df, paste0(save_path, dx_val, "_extcent_diffs.csv"))

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
            n_total_PT = sum(n_total_PT)) %>%
  filter(n_total_CN != 0 & n_total_PT != 0) #if no denominator, skip pheno 
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

fwrite(fisher_final, paste0(save_path, dx_val, "_extcent_fisher.csv"))
