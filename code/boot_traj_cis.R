#Bootstrap CIs for sex-difference trajectories (e.g. cent_plt_pct_mean in review_trajectories.Rmd), for
#each split's curve & the curve averaged across splits:
# - simultaneous bands on centile_M_minus_F_agez (gamlssTools::simultaneous_pct_ci())
# - pointwise CIs on the points of interest (traj_points(): peak ages etc.)

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(tidyr)
library(gamlssTools) #simultaneous_pct_ci()
library(gcplyr) #which_max_gc(), extr_val()

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/"
source(paste0(base, "code/boot_funs.R"))

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
boot_dir_A <- as.character(args[1]) #e.g. cv_sample_A_boot/global_vols_totalTRUE_logAgeTRUE_boot
boot_dir_B <- as.character(args[2]) #matching split B dir
out_file <- as.character(args[3]) #.rds
probs <- c(0.025, 0.5, 0.975)
min_age_days <- 0 #setting low for now
curve_cols <- c("centile_M_minus_F_agez", "deriv_M_minus_F_agez",
                "median_centile_Male", "median_centile_Female",
                "deriv_median_centile_Male", "deriv_median_centile_Female")

phenos <- intersect(basename(list.dirs(boot_dir_A, recursive = FALSE)),
                    basename(list.dirs(boot_dir_B, recursive = FALSE)))
print(paste(length(phenos), "phenos in both splits:", paste(phenos, collapse = ", ")))

read_split <- function(pheno_dir, split) {
  files <- list.files(pheno_dir, pattern = "_uniform_sexdiffs\\.csv$", full.names = TRUE)
  rbindlist(lapply(files, fread, select = c("boot_id", "logAge_days", curve_cols)))[, cv_sample := split]
}

#simultaneous band on one set of curves (one row per boot_id x age_days)
sim_band <- function(curves) {
  #redefine as dataframe
  df <- data.frame(boot = curves$boot_id, age_days = curves$age_days, value = curves$centile_M_minus_F_agez)
  #get CIs
  simultaneous_pct_ci(df, "age_days", "value", interval = probs[3] - probs[1]) %>%
    rename(sim_lower = lower, sim_upper = upper)
}

band_list <- list()
points_list <- list()
n_list <- list()

for (pheno in phenos) {
  print(paste("processing", pheno))
  dt <- rbind(read_split(file.path(boot_dir_A, pheno), "A"),
              read_split(file.path(boot_dir_B, pheno), "B"))
  dt[, age_days := 10^logAge_days] #= gamlssTools::un_log(), as in the Rmd
  dt <- dt[age_days >= min_age_days]

  #each split's curves, plus the split-averaged curves of bootstrap pairs
  #(ages covered by one split only keep that split's value, as in glob_mean)
  pair_ids <- intersect(dt[cv_sample == "A", unique(boot_id)], dt[cv_sample == "B", unique(boot_id)])
  curves <- rbind(
    dt[, c("cv_sample", "boot_id", "age_days", curve_cols), with = FALSE],
    dt[boot_id %in% pair_ids, lapply(.SD, mean), by = .(boot_id, age_days), .SDcols = curve_cols][, cv_sample := "mean"],
    use.names = TRUE)
  setorder(curves, cv_sample, boot_id, age_days)
  curves[, pheno := pheno]

  n_list[[pheno]] <- curves[, .(n_curves = uniqueN(boot_id)), by = .(pheno, cv_sample)]

  #SIMULTANEOUS BANDS on each curve set (A, B, and mean)
  set.seed(1)
  band_list[[pheno]] <- rbindlist(lapply(split(curves, by = "cv_sample"), function(cs) {
    cbind(pheno = pheno, cv_sample = cs$cv_sample[1], sim_band(cs))
  }))

  #POINTS OF INTEREST per bootstrap (A, B) or bootstrap pair (mean)
  points_list[[pheno]] <- traj_points(as.data.frame(curves), by = c("pheno", "cv_sample", "boot_id"))
}

n_curves <- rbindlist(n_list)
bands <- rbindlist(band_list)
boot_points <- bind_rows(points_list)

#POINTWISE CIs on single-valued points (min_sex_diff can have several crossings - see crossing_cis() in the Rmd)
point_cis <- boot_points %>%
  filter(metric != "min_sex_diff") %>%
  summarise(n_boot = sum(!is.na(age)),
            age_lower = quantile(age, probs[1], na.rm = TRUE),
            age_median = quantile(age, probs[2], na.rm = TRUE),
            age_upper = quantile(age, probs[3], na.rm = TRUE),
            value_lower = quantile(value, probs[1], na.rm = TRUE),
            value_median = quantile(value, probs[2], na.rm = TRUE),
            value_upper = quantile(value, probs[3], na.rm = TRUE),
            .by = c(pheno, cv_sample, metric))

print(n_curves)
#n_boot_sim = ceiling(95% of n_curves) by construction; n_tied_* > 1 means the band depends on which tied curves the seed kept
print("simultaneous bands: curves in band & curves tied at the upper/lower cutoffs")
print(bands[, .(n_boot_sim = n_boot_sim[1], n_tied_upper = n_tied_upper[1], n_tied_lower = n_tied_lower[1]),
            by = .(pheno, cv_sample)])

#unscaled values - rescale bands$sim_lower/sim_upper & point_cis$value_* to % of max in the Rmd
print(paste("saving to", out_file))
saveRDS(list(bands = bands, point_cis = point_cis, boot_points = boot_points, n_curves = n_curves,
             probs = probs, min_age_days = min_age_days),
        file = out_file)

print("DONE")
