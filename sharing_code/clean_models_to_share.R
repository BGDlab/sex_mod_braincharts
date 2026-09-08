library(tidyr)
library(data.table)
library(dplyr)
library(purrr)
library(gamlss)
library(gamlssTools)

args <- commandArgs(trailingOnly = TRUE)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "")
mod <- readRDS(args[2])
pheno_cat <- as.character(args[3])
total   <- as.character(args[4])
stopifnot(total %in% c("TRUE", "FALSE"))
split   <- as.character(args[5])
save_path <- as.character(args[6]) #save cleaned models
fig_path <- as.character(args[7]) #save testing outputs

pheno <- get_y(mod) %>% as.character()

#clean model & save
clean_mod <- sanitize_gamlss(mod, grid_n=5000)
clean_filename <- paste0(save_path, "/", pheno, "_split", split, "_total", total, "_sharing_model.rds")
saveRDS(clean_mod, clean_filename)

#compare
z_comp <- compare_scores(mod, clean_mod, data=df, fit_data1=df)
#save z-score comps
z_df <- data.frame(
  phenotype = pheno,
  split = split,
  total = total,
  z_max     = z_comp$z_diffs[["max"]],
  z_mean    = z_comp$z_diffs[["mean"]],
  n_ties    = z_comp$n_tied,
  stringsAsFactors = FALSE
)
z_path <- file.path(fig_path, paste0(".", pheno, "_", split,"_total",total, "_zdiffs.rds"))
saveRDS(z_df, z_path, compress = FALSE)

#plot overlay centile fans
sim_list <- sim_grid(df, "logAge_days", "sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")
plt <- compare_centile_fans(mod,
                            clean_mod,
                            df=df,
                            x_var="logAge_days",
                            facet_var="sexMale",
                            sim_grid_list=sim_list,
                            desiredCentiles = c(0.01, 0.05, 0.5, 0.95, 0.99),
                            x_axis="log_lifespan_fetal") +
  facet_wrap(~sexMale, labeller=labeller(sexMale = c('0'="Female", '1'="Male"))) +
  labs(subtitle=paste("Sample", x$meta$cv_sample, ", total-size-corrected =", x$meta$total, ", ages 2+ only =", x$meta$age2plus))
plt_path <- file.path(fig_path, paste0(".", pheno, "_", split,"_total",total, "_centfan.png"))
ggsave(plt_path, plt)