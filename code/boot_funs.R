#helper functions for bootstrapping test models (fit_boot_mods.R, combine_boot_mods.R)

library(dplyr)

#get_study
#study_site is "<study>_<site>" (see scoring_new_data/data_dictionary.csv) -> keep text before first "_"
#NOTE: assumes study names contain no "_"; check table(get_study(df$study_site)) before running
get_study <- function(study_site) {
  factor(sub("_.*$", "", as.character(study_site)))
}

#draw_boot_idx
#stratified resample with replacement: each stratum keeps its original n, so nrow(boot df) == nrow(df)
#seeded by seed + boot_id so any single bootstrap can be redrawn exactly
draw_boot_idx <- function(df, boot_id, seed, strata = c("study", "study_site")) {
  strata <- match.arg(strata)

  if (strata == "study") {
    grp <- interaction(get_study(df$study_site), df$sexMale, drop = TRUE)
  } else {
    grp <- interaction(df$study_site, df$sexMale, drop = TRUE)
  }

  set.seed(seed + boot_id)
  rows_by_grp <- split(seq_len(nrow(df)), grp)

  #index with sample.int - sample(r) misbehaves when a stratum has a single row
  idx <- unlist(lapply(rows_by_grp, function(r) r[sample.int(length(r), length(r), replace = TRUE)]),
                use.names = FALSE)

  stopifnot(length(idx) == nrow(df))
  return(idx)
}

#subset_start_fv
#gamlss start.from copies *.fv position by position -> reorder test model's fitted values to match boot rows
subset_start_fv <- function(mod, idx) {
  for (p in mod$parameters) {
    fv_name <- paste0(p, ".fv")
    mod[[fv_name]] <- mod[[fv_name]][idx]
  }
  return(mod)
}

#simultaneous_ci
#simultaneous percentile CIs over age (Gao, Konietschke & Li 2021, doi:10.3390/sym13071212)
#rewrite of gamlssTools::simultaneous_pct_ci()
#  dt: long data.table, one row per bootstrap x age, with columns boot_id, value & age_var
#returns per age: band bounds, n_boot_sim (curves in band), & n_tied_upper/n_tied_lower (curves tied at each
#cutoff - when > 1, which tied curves are kept depends on the seed)
simultaneous_ci <- function(dt, age_var, interval = 0.95) {
  alpha <- 1 - interval
  dt <- dt[!is.na(value)]
  n_boot <- uniqueN(dt$boot_id)
  n_phi <- ceiling(n_boot * (1 - alpha / 2)) #curves kept after the upper trim
  n_final <- ceiling(n_boot * (1 - alpha)) #curves kept in the band

  #upper trim: keep the n_phi bootstraps with the lowest highest-rank across ages (ties at the cutoff split at random)
  dt[, rank := frank(value, ties.method = "random"), by = c(age_var)]
  max_rank <- dt[, .(rank_stat = max(rank)), by = boot_id][order(rank_stat, runif(.N))]
  n_tied_upper <- max_rank[rank_stat == rank_stat[n_phi], .N]
  phi <- dt[boot_id %in% max_rank$boot_id[seq_len(n_phi)]]

  #lower trim: re-rank within phi & keep the n_final bootstraps with the highest lowest-rank (ties split at random)
  phi[, rank := frank(value, ties.method = "random"), by = c(age_var)]
  min_rank <- phi[, .(rank_stat = min(rank)), by = boot_id][order(-rank_stat, runif(.N))]
  n_tied_lower <- min_rank[rank_stat == rank_stat[n_final], .N]
  final <- phi[boot_id %in% min_rank$boot_id[seq_len(n_final)]]

  #band = envelope of the remaining curves
  band <- final[, .(sim_lower = min(value), sim_upper = max(value), n_boot_sim = uniqueN(boot_id)), by = c(age_var)]
  band[, `:=`(n_tied_upper = n_tied_upper, n_tied_lower = n_tied_lower)]
}
