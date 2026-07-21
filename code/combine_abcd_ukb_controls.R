# Combine per-phenotype control centiles (from centile_calc.R outputs) into a
# single wide CSV, keeping only ABCD + UKB controls and only _full-model centiles.
#
# centile_calc.R writes controls and patients to separate files:
#   <save_path>/cent_csvs/<pheno>_CN_<dx>_cent.csv   (controls)
#   <save_path>/cent_csvs/<pheno>_PT_<dx>_cent.csv   (patients)
# so we only need the _CN_ files. Each file is one phenotype (centile columns
# named <pheno>_*), and each phenotype has up to 6 _CN_ files (one per dx). A
# control subject's centile for a phenotype is identical across dx comparisons,
# so we union the dx files per phenotype and dedup by INDEX.ID, then join wide.

library(data.table); library(purrr)

base_dir  <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts"
studies   <- c("ABCD", "UKB")   # <-- match your exact study labels
study_col <- "study"            # <-- CONFIRM: likely "study" or "study_site"

# controls only, main analyses only (drop weighted/age2plus sensitivity runs)
cn_files <- Sys.glob(file.path(base_dir, "cv_sample_*_test",
                               "*totalTRUE*logAgeTRUE_pbmods",
                               "cent_csvs", "*_CN_*_cent.csv"))
cn_files <- cn_files[!grepl("weighted|age2plus", cn_files)]

# keep INDEX.ID + this phenotype's _full centile columns, ABCD/UKB rows only
read_cn <- function(f) {
  pheno <- sub("_CN_.*", "", basename(f))
  dt <- fread(f)[get(study_col) %in% studies]
  cent_cols <- grep(paste0("^", pheno, "_.*_full$"), names(dt), value = TRUE)
  dt[, c("INDEX.ID", cent_cols), with = FALSE]
}

# per phenotype: union controls across dx comparisons, dedup by subject
per_pheno <- split(cn_files, sub("_CN_.*", "", basename(cn_files))) |>
  map(\(fs) unique(rbindlist(map(fs, read_cn)), by = "INDEX.ID"))

# wide join: one row per control, all phenotype centiles
combined <- reduce(per_pheno, \(x, y) merge(x, y, by = "INDEX.ID", all = TRUE))

# reattach shared covariates once
covars <- unique(rbindlist(lapply(cn_files, \(f) {
  dt <- fread(f)[get(study_col) %in% studies]
  keep <- intersect(c("INDEX.ID", "age_days", "sex", study_col, "dx", "dx_recode"), names(dt))
  dt[, ..keep]
})), by = "INDEX.ID")

combined <- merge(covars, combined, by = "INDEX.ID", all.y = TRUE)
fwrite(combined, file.path(base_dir, "abcd_ukb_controls_all_phenotypes_cent.csv"))
