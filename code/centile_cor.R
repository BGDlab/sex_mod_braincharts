library(tidyr)
library(data.table)
library(dplyr)
library(purrr)
library(EnvStats)

#set paths
base_path <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
save_path <- paste0(base_path, "cent_cor/")

#define the grid of conditions to run
dx_levels    <- c("ADHD", "ALZ", "ASD", "GAD", "MDD", "SCZ")
total_levels <- c("TRUE", "FALSE")

#get pheno list
lists <- list.files(paste0(base_path, "pheno_lists"), pattern = "\\.rds$", full.names = TRUE)
pheno_list <- do.call(c, lapply(lists, readRDS))

fread_filt <- function(f, string, split){
  fread(f) %>%
    select(INDEX.ID, sexMale, dx_recode, logAge_days, matches(string)) %>%
    mutate(sex = ifelse(sexMale==0, "F", "M"),
           split = split)
}

#compute correlations for one (dx_val, total) combination
run_combo <- function(dx_val, total) {
  csv_path <- paste0(base_path, "cv_sample_?_test/*total", total, "*")
  
  results_list <- list()
  for (pheno in pheno_list) {
    f_list <- Sys.glob(paste0(csv_path, "/cent_csvs/", pheno, "_PT_", dx_val, "_cent.csv"))
    if (length(f_list) != 2) {
      warning(paste(length(f_list), "file(s) found for", pheno, "dx", dx_val,
                    "total", total, "- skipping"))
      next
    }
    names(f_list) <- sub(".*cv_sample_(.).*", "\\1", f_list)
    pheno_df <- rbindlist(
      Map(fread_filt, f_list, "full", names(f_list)),
      fill = TRUE
    )
    
    cent_col <- grep("_centile_full$",   names(pheno_df), value = TRUE)
    std_col  <- grep("_std_score_full$", names(pheno_df), value = TRUE)
    
    run_cor <- function(col, measure) {
      wide <- dcast(pheno_df, INDEX.ID ~ split, value.var = col)
      ct <- cor.test(wide$A, wide$B)
      data.table(
        dx      = dx_val,
        total   = total,
        pheno   = pheno,
        measure = measure,
        column  = col,
        r       = ct$estimate,
        p       = ct$p.value,
        ci_low  = ct$conf.int[1],
        ci_high = ct$conf.int[2],
        n       = sum(complete.cases(wide$A, wide$B))
      )
    }
    
    results_list[[pheno]] <- rbind(
      run_cor(cent_col, "centile"),
      run_cor(std_col,  "std_score")
    )
  }
  
  if (length(results_list) == 0) return(NULL)
  results <- rbindlist(results_list)
}

#run over the full grid
grid <- expand.grid(dx = dx_levels, total = total_levels, stringsAsFactors = FALSE)
all_results <- rbindlist(
  Map(run_combo, grid$dx, grid$total),
  fill = TRUE
)

#combined output
fwrite(all_results, paste0(save_path, "split_correlations_ALL.csv"))

p <- all_results %>%
  group_by(total, dx, measure) %>%
  summarise(mean_r=mean(r)) %>%
  ggplot() +
  geom_col(aes(x=measure, y=mean_r, fill=dx), position="dodge") +
  facet_wrap(~total)

p <- all_results %>%
  group_by(total, dx, measure) %>%
  ggplot() +
  geom_boxplot(aes(x=measure, y=r, fill=dx), position="dodge") +
  facet_wrap(~total)
