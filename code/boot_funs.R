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


#Functions for the points of interest on sex-difference trajectories 
#requires dplyr, tidyr, gcplyr

#traj_points
#same metrics as glob_mean_sums + min_sex_diff_rows in review_trajectories.Rmd, one set per group in `by`
#  curves: one row per group x age, sorted by age within group, with age_days, centile_M_minus_F_agez,
#          deriv_M_minus_F_agez, median_centile_Male/Female & deriv_median_centile_Male/Female
#  by: grouping columns, e.g. "pheno" for the actual curves, c("pheno", "cv_sample", "boot_id") for bootstraps
#returns long df: by columns, metric, value, age, crossing (TRUE for min_sex_diff rows at a zero crossing,
#FALSE for the closest-to-zero fallback when a curve never crosses 0, NA for other metrics)
traj_points <- function(curves, by = "pheno") {
  sums <- curves %>%
    summarize(#age at peak sex diff
      max_sex_diff = extr_val(centile_M_minus_F_agez, which_max_gc(abs(centile_M_minus_F_agez))),
      age_max_sex_diff = extr_val(age_days, which_max_gc(abs(centile_M_minus_F_agez))),
      #age peak change in sex diff
      max_deriv = extr_val(centile_M_minus_F_agez, which_max_gc(abs(deriv_M_minus_F_agez))),
      age_max_deriv = extr_val(age_days, which_max_gc(abs(deriv_M_minus_F_agez))),
      #age at male peak
      max_M = extr_val(centile_M_minus_F_agez, which_max_gc(median_centile_Male)),
      age_max_M = extr_val(age_days, which_max_gc(median_centile_Male)),
      #age at female peak
      max_F = extr_val(centile_M_minus_F_agez, which_max_gc(median_centile_Female)),
      age_max_F = extr_val(age_days, which_max_gc(median_centile_Female)),
      #age at male peak growth
      max_deriv_M = extr_val(centile_M_minus_F_agez, which_max_gc(deriv_median_centile_Male)),
      age_max_deriv_M = extr_val(age_days, which_max_gc(deriv_median_centile_Male)),
      #age at female peak growth
      max_deriv_F = extr_val(centile_M_minus_F_agez, which_max_gc(deriv_median_centile_Female)),
      age_max_deriv_F = extr_val(age_days, which_max_gc(deriv_median_centile_Female)),
      .by = all_of(by))
  
  #pivot and rejoin
  values_long <- sums %>%
    select(!starts_with("age_")) %>%
    pivot_longer(cols = -all_of(by), names_to = "metric", values_to = "value")
  
  ages_long <- sums %>%
    select(all_of(by), starts_with("age_")) %>%
    pivot_longer(cols = starts_with("age_"), names_to = "metric", values_to = "age") %>%
    mutate(metric = gsub("^age_", "", metric))
  
  points_long <- left_join(values_long, ages_long, by = c(by, "metric"))
  
  #age at min sex-diff: all zero crossings, taking whichever flanking row has smaller |value|;
  #curves with no crossing fall back to the single closest-to-zero row
  min_sex_diff_rows <- curves %>%
    arrange(across(all_of(by)), age_days) %>%
    group_by(across(all_of(by))) %>%
    mutate(
      sign_now  = sign(centile_M_minus_F_agez),
      sign_next = sign(lead(centile_M_minus_F_agez)),
      crosses   = !is.na(sign_next) & sign_now != 0 & sign_next != 0 & sign_now != sign_next,
      # mark the closer-to-zero row of each adjacent crossing pair
      pick_here = crosses & abs(centile_M_minus_F_agez) <= abs(lead(centile_M_minus_F_agez)),
      pick_next = lag(crosses & abs(centile_M_minus_F_agez) >  abs(lead(centile_M_minus_F_agez)),
                      default = FALSE),
      keep      = pick_here | pick_next,
      any_cross = any(keep, na.rm = TRUE),
      # fallback: if no crossing, keep the single global min |value|
      keep      = if (any_cross[1]) keep
      else abs(centile_M_minus_F_agez) == min(abs(centile_M_minus_F_agez), na.rm = TRUE)
    ) %>%
    filter(keep) %>%
    ungroup() %>%
    transmute(
      across(all_of(by)),
      metric   = "min_sex_diff",
      value    = centile_M_minus_F_agez,
      age      = age_days,
      crossing = any_cross
    )
  
  bind_rows(points_long, min_sex_diff_rows)
}

#crossing_cis
#CIs on the ages of the actual curve's min_sex_diff points. Bootstraps can cross 0 a different number of times,
#so each bootstrap crossing is matched to the nearest actual crossing (on log age), keeping a bootstrap's closest
#crossing per actual crossing. pct_boot = % of bootstrap curves with a crossing matched to it - low values mean a
#crossing many bootstraps don't reproduce.
#If the actual curve never crosses 0, its fallback (closest-to-zero) point gets a CI from each bootstrap's
#closest-to-zero min_sex_diff row, & pct_boot = % of bootstrap curves that also never cross.
#  boot_points: boot_traj_cis.R's $boot_points (traj_points(..., by = c(by, "boot_id")))
#  real_points: traj_points(..., by = by) on the actual curves, e.g. glob_mean_sums_long w/ cv_sample = "mean"
#  n_curves: df of `by` columns & n_curves (boot_traj_cis.R's $n_curves)
#  by: columns identifying a curve, e.g. c("pheno", "cv_sample") for boot_traj_cis.R's split & averaged curves
crossing_cis <- function(boot_points, real_points, n_curves, by = c("pheno", "cv_sample"),
                         probs = c(0.025, 0.5, 0.975)) {
  real <- real_points %>%
    filter(metric == "min_sex_diff") %>%
    select(all_of(by), real_age = age, real_crossing = crossing) %>%
    mutate(crossing_id = row_number(), .by = all_of(by))
  boot <- boot_points %>%
    filter(metric == "min_sex_diff") %>%
    semi_join(real, by = by) %>%
    select(all_of(by), boot_id, age, value, crossing)
  
  #actual curve crosses 0: match bootstrap crossings to the nearest actual crossing
  matched_cross <- boot %>%
    filter(crossing) %>%
    inner_join(filter(real, real_crossing), by = by, relationship = "many-to-many") %>%
    mutate(dist = abs(log10(age) - log10(real_age))) %>%
    slice_min(dist, n = 1, with_ties = FALSE, by = all_of(c(by, "boot_id", "age"))) %>%
    slice_min(dist, n = 1, with_ties = FALSE, by = all_of(c(by, "boot_id", "crossing_id"))) %>%
    select(all_of(by), boot_id, crossing_id, real_age, age)
  
  #actual curve never crosses 0: each bootstrap's closest-to-zero row
  matched_fallback <- boot %>%
    inner_join(filter(real, !real_crossing), by = by) %>%
    slice_min(abs(value), n = 1, with_ties = FALSE, by = all_of(c(by, "boot_id"))) %>%
    select(all_of(by), boot_id, crossing_id, real_age, age, crossing)
  
  support_fallback <- matched_fallback %>%
    summarise(n_support = sum(!crossing), .by = all_of(c(by, "crossing_id")))
  
  bind_rows(matched_cross, select(matched_fallback, !crossing)) %>%
    summarise(n_boot = n(),
              age_lower = quantile(age, probs[1]),
              age_median = quantile(age, probs[2]),
              age_upper = quantile(age, probs[3]),
              .by = all_of(c(by, "crossing_id", "real_age"))) %>%
    right_join(select(real, all_of(by), crossing_id, real_age, real_crossing),
               by = c(by, "crossing_id", "real_age")) %>%
    left_join(support_fallback, by = c(by, "crossing_id")) %>%
    left_join(n_curves, by = by) %>%
    mutate(n_boot = coalesce(n_boot, 0L),
           pct_boot = if_else(real_crossing, n_boot, n_support) / n_curves * 100) %>%
    select(!n_support) %>%
    arrange(across(all_of(by)), real_age)
}