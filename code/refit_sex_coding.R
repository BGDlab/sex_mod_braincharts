#Refit a training BestMod (sexMale: female=0, male=1) on the same data with sexFemale (male=0, female=1)
#and sexEffect (female=-1, male=1) coding, and with pb(logAge_days) + pvc(logAge_days, by = sexEffect) in place of
#the sex x age smooth, then plot centile fans + sex-specific worm plots for all 4 codings.
#Outputs go to <model dir>/sex_coding/<model name>/ and are compiled by compile_sex_coding_pdf.R
set.seed(99999)

#LOAD PACKAGES
library(data.table)
library(ggplot2)
library(dplyr)
library(gamlss)
library(gamlssTools)
library(ggnewscale)

base <- Sys.getenv("SEX_MOD_BASE", "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/")
source(paste0(base, "code/gamlss_fit_funs.R"))

#centile fans predict through gamlssTools, which rebuilds pvc() terms from the fit: gamlss's own predict(newdata =)
#centres pvc's numeric `by` on old + new data together, shifting every prediction
if (!exists(".pvc_spec", envir = asNamespace("gamlssTools"), inherits = FALSE))
  stop("gamlssTools ", packageVersion("gamlssTools"), " has no pvc() support; ",
       "install the dev branch (commit c080639 or later) so the pvc coding plots correctly")

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
mod_file <- normalizePath(args[1]) #training BestMod (.rds)
csv <- if (length(args) >= 2 && !(args[2] %in% c("", "NA"))) args[2] else NA #optional: training csv, else found from mod_file path

og_mod <- readRDS(mod_file)
pheno <- as.character(og_mod$mu.terms[[2]])
stopifnot("BestMod must use sexMale_x_logAge" = get_sex_var(deparse1(og_mod$mu.formula)) == "sexMale" &&
            grepl("sexMale_x_logAge", deparse1(og_mod$mu.formula), fixed = TRUE))

#model dir: <root>/cv_sample_<split>_train/<pheno_cat>_total<T>_logAge<L>_pbmods/model_objs/<pheno>_<m_name>_BestMod.rds
save_dir <- dirname(dirname(mod_file))
mod_set <- basename(save_dir)
total_tf <- regmatches(mod_set, regexpr("(?<=_total)(TRUE|FALSE)", mod_set, perl = TRUE))
log_age <- regmatches(mod_set, regexpr("(?<=_logAge)(TRUE|FALSE)", mod_set, perl = TRUE))
sample_name <- sub("_train$", "", basename(dirname(save_dir)))
m_name <- sub("_BestMod\\.rds$", "", basename(mod_file))

#find training csv the same way config_cv_mods.sh does: data/<sample>_dfs/<pheno>_total<T>*logAge<L>.csv
if (is.na(csv)) {
  csv_dir <- file.path(dirname(dirname(save_dir)), "data", paste0(sample_name, "_dfs"))
  hits <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  hits <- hits[startsWith(basename(hits), paste0(pheno, "_total", total_tf)) &
                 endsWith(basename(hits), paste0("logAge", log_age, ".csv"))]
  if (length(hits) != 1) stop("expected 1 training csv for ", pheno, " in ", csv_dir, ", found ", length(hits))
  csv <- hits
}
print(paste("training data:", csv))

df_og <- droplevels(fread(csv, stringsAsFactors = TRUE, na.strings = c("", '""'))) #same read as fit_cv_mods.R

#make sure this is the data the BestMod was fit on
if (nrow(df_og) != og_mod$N || !isTRUE(all.equal(as.numeric(df_og[[pheno]]), as.numeric(og_mod$y), check.attributes = FALSE))) {
  stop(csv, " does not match the data ", mod_file, " was fit on (", nrow(df_og), " vs ", og_mod$N, " rows)")
}

out_dir <- file.path(save_dir, "sex_coding", m_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

coding_labels <- c(female_ref = "sexMale (F=0, M=1), original",
                   male_ref = "sexFemale (M=0, F=1)",
                   effect = "sexEffect (F=-1, M=1)",
                   pvc = "pb(logAge) + pvc(logAge, by = sexEffect) (F=-1, M=1)")
total_candidates <- c("TBV", "mean.CT", "total.SA")
birth <- log(280, base=10)

stats_df <- data.frame()

for (sex_coding in names(sex_codings)) {
  sex_var <- sex_codings[[sex_coding]]
  sex_age_var <- paste0(sex_var, "_x_logAge")
  sex_age_lab <- sexage_label(sex_coding) #how the sex x age smooth appears in the model's term labels
  print(paste("sex coding:", sex_coding, "->", sex_var))

  #global df: model formulas use k=log(nrow(df)) and data=df
  df <- if (sex_coding == "female_ref") df_og else add_sex_coding(df_og, sex_coding)

  #GET MODEL
  if (sex_coding == "female_ref") {
    model <- og_mod
  } else {
    #(effect & pvc both use sexEffect, so pvc gets its own file name)
    m_file <- file.path(out_dir, paste0(m_name, "_", if (sex_coding == "pvc") "pvc_sexEffect" else sex_var, ".rds"))
    model <- NULL
    if (file.exists(m_file)) {
      print("loading pre-fit model")
      model <- tryCatch(readRDS(m_file), error = function(e) NULL)
      if (!is.null(model) && (!isTRUE(model$converged) || get_sex_coding(deparse1(model$mu.formula)) != sex_coding)) {
        print("pre-fit model not converged or wrong coding - refitting")
        model <- NULL
      }
    }
    if (is.null(model)) {
      print("refitting")
      model <- refit_sex_coding(og_mod, sex_coding)
      if (!is.null(model)) saveRDS(model, file = m_file)
    }
  }

  if (is.null(model)) {
    message("refit failed for ", sex_coding)
    stats_df <- rbind(stats_df, data.frame(coding = sex_coding, sex_var = sex_var, converged = FALSE))
    next
  }
  model$call$data <- "df"
  model$call$family <- model$family[[1]]

  #STATS
  #edf & lambda of the sex x age smooth in mu & sigma (lambda 1e7 = capped, i.e. linear in age)
  smo_stat <- function(what, stat) {
    smo <- model[[paste0(what, ".coefSmo")]]
    i <- grep(sex_age_lab, colnames(model[[paste0(what, ".s")]]), fixed = TRUE)
    if (length(i) == 1) smo[[i]][[stat]] else NA
  }
  smo_edf <- function(what) smo_stat(what, "edf")
  stats_df <- rbind(stats_df, data.frame(
    coding = sex_coding,
    sex_var = sex_var,
    converged = isTRUE(model$converged),
    deviance = deviance(model),
    df_fit = model$df.fit,
    AIC = model$aic,
    BIC = model$sbc,
    mu_sexAge_edf = smo_edf("mu"),
    sigma_sexAge_edf = smo_edf("sigma"),
    mu_sexAge_lambda = smo_stat("mu", "lambda"),
    sigma_sexAge_lambda = smo_stat("sigma", "lambda")
  ))

  #CENTILE FAN PLOT (as in replot_centiles.R)
  pred_list <- setdiff(list_predictors(model), "df")
  special_term <- sex_coding_special_terms(sex_coding)
  total_var <- intersect(total_candidates, pred_list)
  if (length(total_var) > 1) stop("multiple total vars in model: ", paste(total_var, collapse = ", "))
  if (length(total_var) == 1) special_term <- c(special_term, paste0(total_var, "=0")) #hold total at 0 for plotting
  pt_show <- length(total_var) == 0

  print("simulate data for plotting")
  sim_df <- sim_data(df, "logAge_days", factor_var = "sexMale", special_term = special_term)
  resid_terms <- setdiff(pred_list, c(sex_age_var, "logAge_days", sex_var))

  fan <- make_centile_fan(gamlssModel = model,
                          df = df,
                          x_var = "logAge_days",
                          color_var = "sexMale",
                          get_peaks = TRUE,
                          desiredCentiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                          show_points = pt_show,
                          sim_data_list = sim_df,
                          label_centiles = "legend",
                          point_color_manual = c('0' = "#F4A15BFF", '1' = "#8CB3D1FF"),
                          color_manual = c("0" = "#c05f0d", "1" = "#38688D"),
                          remove_point_effect = resid_terms,
                          x_axis = "log_lifespan_fetal") +
    theme_linedraw() +
    scale_color_manual(name = "Sex",
                       values = c('0' = "#c05f0d", '1' = "#38688D"),
                       labels = c("Female", "Male")) +
    guides(fill = "none",
           fill_ggnewscale_1 = "none", colour_ggnewscale_1 = "none") +
    geom_vline(xintercept = birth, color = "gray20", linetype = "dashed") +
    labs(title = paste0(pheno, ": ", coding_labels[[sex_coding]]),
         subtitle = if (length(total_var) == 1) paste(total_var, "held at 0") else NULL,
         x = "Age (years)")
  ggsave(file.path(out_dir, paste0("fan_", sex_coding, ".png")), fan, width = 8, height = 5)

  #SEX-SPECIFIC WORM PLOTS
  res <- resid(model)
  for (s in c("0", "1")) {
    idx <- as.character(df$sexMale) == s
    wp <- wp.taki(xvar = df$logAge_days[idx], resid = res[idx], n.inter = 6)$plot +
      ggtitle(paste0(ifelse(s == "0", "Female", "Male"), " residuals: ", coding_labels[[sex_coding]]))
    ggsave(file.path(out_dir, paste0("worm_", sex_coding, "_", ifelse(s == "0", "F", "M"), ".png")), wp, width = 6, height = 5)
  }

  rm(model)
}

stats_df$pheno <- pheno
stats_df$m_name <- m_name
stats_df$sample <- sample_name
stats_df$total <- total_tf
stats_df$bestmod <- mod_file
stats_df$csv <- csv
stats_df$n <- nrow(df_og)
fwrite(stats_df, file.path(out_dir, "sex_coding_stats.csv"))
print(stats_df)

print("DONE")
