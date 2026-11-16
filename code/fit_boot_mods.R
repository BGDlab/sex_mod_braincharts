#Fit one bootstrap replicate of a test model:
#stratified resample of the test df, refit with the training model's fixed lambdas, warm start from the test model

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc/"
source(paste0(base, "code/gamlss_fit_funs.R"))
source(paste0(base, "code/boot_funs.R"))
source(paste0(base, "code/sex_diff_funs.R"))
options(warn = 1)

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df_full <- droplevels(fread(args[1], stringsAsFactors = TRUE, na.strings = c("", '""'))) #same read as fit_cv_mods_test.R
base_mod <- readRDS(args[2]) #training BestMod -> fixed lambdas
test_path <- as.character(args[3]) #test full_mod -> warm start
save_path <- as.character(args[4])
sim_df_uniform <- readRDS(args[5])
boot_id <- as.integer(args[6])
seed <- as.integer(args[7])
strata <- as.character(args[8])
n.cyc <- 1800

#SEX CODING: pb(logAge_days) + pvc(logAge_days, by = sexEffect) (see gamlss_fit_funs.R) - same checks as fit_cv_mods_test.R
stopifnot("training model has no pvc sex x age term" = grepl("pvc(logAge_days", deparse1(base_mod$mu.formula), fixed = TRUE))
stopifnot("df has no sexEffect column (rerun lbcc_eda.Rmd)" = SEX_VAR %in% names(df_full),
          "sim_df_uniform's sexEffect doesn't match sexMale (rebuild v3_full_sim.rds with special_term = SEX_SPECIAL_TERM)" =
            all(vapply(sim_df_uniform, function(s) SEX_VAR %in% names(s) && all(s[[SEX_VAR]] == 2 * s$sexMale - 1), logical(1))))

pheno <- base_mod$mu.terms[[2]] %>% as.character()
prefix <- paste0(save_path, "/boot_", sprintf("%03d", boot_id))
file_diffs <- paste0(prefix, "_uniform_sexdiffs.csv")
file_info <- paste0(prefix, "_info.csv")
file_coefs <- paste0(prefix, "_coefs.csv")

#skip if already done
if (file.exists(file_diffs) & file.exists(file_info)) {
  print(paste("already fit:", file_diffs))
  quit(save = "no", status = 0)
}

#RESAMPLE
#df & test_mod must be globals named exactly this: gamlss_lambda_rep() eval(parse())'s "data = df" & "start.from = test_mod"
idx <- draw_boot_idx(df_full, boot_id, seed, strata)
df <- droplevels(df_full[idx, ])

print("strata (original df):")
print(table(df_full[[strata]], df_full$sexMale))

#WARM START
test_mod <- readRDS(test_path)
stopifnot(length(test_mod$mu.fv) == nrow(df_full)) #fv must be in df_full row order for idx subsetting to be valid
test_mod <- subset_start_fv(test_mod, idx)
stopifnot(length(test_mod$mu.fv) == nrow(df))

#FIT
base_mod$call$data <- "df"
base_mod$call$family <- base_mod$family[[1]]

age_var <- if ("logAge_days" %in% names(df)) "logAge_days" else "age_days"

info_df <- data.frame(
  "pheno" = pheno,
  "boot_id" = boot_id,
  "seed" = seed + boot_id,
  "strata" = strata,
  "n" = nrow(df),
  "n_unique" = length(unique(idx)),
  "idx_hash" = rlang::hash(idx), #fingerprint of the resample - compare across reruns/scripts to confirm identical draws
  "age_min" = min(df[[age_var]]),
  "age_max" = max(df[[age_var]]),
  "n_sites" = nlevels(df$study_site),
  "n_sites_full" = nlevels(df_full$study_site)
)

t0 <- Sys.time()
model <- gamlss_lambda_rep(base_mod, null_mod = "false", start.from = "test_mod", n.cyc = n.cyc)
info_df$fit_mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

if (is.null(model)) {
  info_df$converged <- FALSE
  info_df$iter <- NA
  fwrite(info_df, file = file_info)
  stop("model fitting failed")
}

info_df$converged <- model$converged
info_df$iter <- model$iter

#lambdas are fixed by construction - check that every smoother kept them (no leftover GAIC terms)
for (p in c("mu", "sigma", "nu")) {
  ref <- base_mod[[paste0(p, ".lambda")]]
  if (!is.null(ref)) stopifnot(max(abs(ref - model[[paste0(p, ".lambda")]])) < 0.01)
}

model$call$data <- "df"
model$call$family <- model$family[[1]]

#COEFS
coef_df <- rbindlist(lapply(model$parameters, function(p) {
  cf <- coef(model, what = p)
  data.table(boot_id = boot_id, moment = p, term = names(cf), estimate = unname(cf))
}))

#SEX DIFFS on uniform sim grid
result_df <- process_sex_diffs(sim_data_list = sim_df_uniform, model = model, df = df, age_var = age_var,
                               age_range = range(df_full[[age_var]])) %>%
  mutate(boot_id = boot_id)

#SAVE - info last so its presence marks a finished job
print(paste("saving to", prefix))
fwrite(coef_df, file = file_coefs)
fwrite(result_df, file = file_diffs)
fwrite(info_df, file = file_info)

print("SUCCESS")
