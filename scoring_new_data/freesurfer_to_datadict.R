#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# freesurfer_to_datadict.R
#
# Reshape FreeSurfer output into the wide, one-row-per-scan CSV for scoring described by
# data_dictionary.csv:
#   68 Desikan-Killiany regions x {SA, GM, CT}  (aparc.stats)
#   30 subcortical volumes  (aseg.stats)
#   global tissue volumes GMV / sGMV / WMV / CSF / CBV / TBV, mean.CT, total.SA
#   fs_version_SA / _CT / _GM, study_site, sexMale, logAge_days
#
# WARNING: Demographics matching assumes cross-sectional data (1 unique demographics row per participant).
# No QC filtering
#
# Two input modes:
#
# 1. recon-all subject directories (default) -- reads
#      <subj>/stats/{lh.aparc.stats, rh.aparc.stats, aseg.stats}
#
#      Rscript scoring_new_data/freesurfer_to_datadict.R \
#        --subjects-dir /path/to/SUBJECTS_DIR \
#        --demographics demo.csv --demo-id-col participant \
#        --age-col age_days --sex-col sex \
#        --study ABCD --site-col site \
#        --fs-version FS7_T1 \
#        --out abcd_datadict.csv
#
# 2. group tables from aparcstats2table / asegstats2table
#
#      Rscript scoring_new_data/freesurfer_to_datadict.R \
#        --lh-area lh.area.tsv --rh-area rh.area.tsv \
#        --lh-volume lh.vol.tsv --rh-volume rh.vol.tsv \
#        --lh-thickness lh.thick.tsv --rh-thickness rh.thick.tsv \
#        --aseg aseg.tsv \
#        --demographics demo.csv --demo-id-col participant \
#        --age-col age_years --age-units years --sex-col sex \
#        --fs-version FS6_T1 \
#        --out abcd_datadict.csv
#
# Demographics are joined on --demo-id-col, matched against the FreeSurfer subject
# ID. Everything else is computed from the stats.
#
# Mappings:
#   SA = aparc SurfArea, GM = aparc GrayVol, CT = aparc ThickAvg
#   GMV  = CortexVol, sGMV = SubCortGrayVol,
#   WMV  = CerebralWhiteMatterVol (CorticalWhiteMatterVol in FS 5.3)
#   CSF  = BrainSegVol - BrainSegVolNotVent
#   CBV  = L/R cerebellum cortex + L/R cerebellum white matter
#   TBV  = WMV + sGMV + GMV + CBV
#   mean.CT / total.SA = unweighted mean / sum over the 68 DK regions; NA if any
#     region is missing
#   SUBC.*.Thalamus.Proper accepts the FS7 "Left-Thalamus" / "Right-Thalamus" name
#   logAge_days = log10(age_days + 280) unless --age-type postconception
#     (280 is overridable with a per-subject --ga-col)
#   sexMale = 1/0 via --male-values / --female-values; anything else is NA
#
# Options (written as "--key value" or "--key=value"):
#
#   input -- recon-all directories
#     --subjects-dir DIR     SUBJECTS_DIR holding <subj>/stats/*.stats
#     --subjects-list FILE   subject IDs to include, one per line, each resolved
#                            as <subjects-dir>/<id>. If NULL (default) every directory
#                            under --subjects-dir that has stats/aseg.stats.
#                            List can be used to convert a subset, to fix the row
#                            order, to split a cohort across array jobs, or to
#                            make failed subjects show up as an all-NA row
#                            instead of silently vanishing. Only read alongside
#                            --subjects-dir; ignored in table mode.
#
#   input -- aparcstats2table / asegstats2table group tables (use instead of --subjects-dir)
#     --lh-area FILE         aparcstats2table --hemi lh --meas area
#     --rh-area FILE           "                   rh
#     --lh-volume FILE       aparcstats2table --hemi lh --meas volume
#     --rh-volume FILE         "                   rh
#     --lh-thickness FILE    aparcstats2table --hemi lh --meas thickness
#     --rh-thickness FILE      "                   rh
#     --aseg FILE            asegstats2table --meas volume
#
#   demographics (a CSV joined to the FreeSurfer subject IDs; necessary for
#   demographic covariates)
#     --demographics FILE    the CSV
#     --demo-id-col NAME          its column matched against the subject ID;
#                            defaults to "participant"
#     --age-col NAME         age at scan -> logAge_days
#     --age-units UNIT       days (default) | weeks | months | years
#     --age-type TYPE        birth (default; 280 or --ga-col value is added) |
#                            postconception (used as-is)
#     --ga-col NAME          per-subject gestational age at birth, overriding
#                            280 days where it is not NA
#     --ga-units UNIT        units of --ga-col: days (default) | weeks |
#                            months | years
#     --sex-col NAME         sex -> sexMale
#     --male-values CSV      values recoded to 1, default
#                            "M,m,Male,male,MALE"
#     --female-values CSV    values recoded to 0, default
#                            "F,f,Female,female,FEMALE"
#                            (matching is exact after trimming whitespace;
#                            anything in neither list becomes NA and is warned
#                            about once per distinct value)
#
#   labels written into every row
#     --study NAME           study label; with --site-col it prefixes the site
#     --site-col NAME        per-subject site column in the demographics CSV;
#                            study_site becomes "<study>_<site>"
#     --fs-version V         fills fs_version_SA, _CT and _GM; see 'data_dictionary.csv'
#     --fs-version-sa V      override just fs_version_SA
#     --fs-version-ct V      override just fs_version_CT
#     --fs-version-gm V      override just fs_version_GM
#                            (V should be one of FS_VERSIONS included in training data
#                             (listed below); anything else is written through with a 
#                            warning. See 'data_dictionary.csv' for more info)
#
#   output
#     --out FILE             required; parent directories are created
#     --id-name NAME         name of the leading subject-ID column, default
#                            "participant"
#
#   --help / -h              print the usage summary and the option names
# ---------------------------------------------------------------------------

# ---- column definitions ----------------

DK_REGIONS <- c(
  "bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", "cuneus",
  "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal",
  "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", "lingual",
  "medialorbitofrontal", "middletemporal", "parahippocampal", "paracentral",
  "parsopercularis", "parsorbitalis", "parstriangularis", "pericalcarine",
  "postcentral", "posteriorcingulate", "precentral", "precuneus",
  "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal",
  "superiorparietal", "superiortemporal", "supramarginal", "frontalpole",
  "temporalpole", "transversetemporal", "insula"
)

# dictionary names (minus the SUBC. prefix); FreeSurfer StructName is the same
# string with "." -> "-", except for thalamus (below)
SUBC_STRUCTS <- c(
  "Left.Lateral.Ventricle", "Left.Cerebellum.White.Matter",
  "Left.Cerebellum.Cortex", "Left.Thalamus.Proper", "Left.Caudate",
  "Left.Putamen", "Left.Pallidum", "3rd.Ventricle", "4th.Ventricle",
  "Brain.Stem", "Left.Hippocampus", "Left.Amygdala", "Left.Accumbens.area",
  "Left.VentralDC", "Right.Cerebellum.White.Matter", "Right.Cerebellum.Cortex",
  "Right.Thalamus.Proper", "Right.Caudate", "Right.Putamen", "Right.Pallidum",
  "Right.Hippocampus", "Right.Amygdala", "Right.Accumbens.area",
  "Right.VentralDC", "Right.Lateral.Ventricle", "CC_Posterior",
  "CC_Mid_Posterior", "CC_Central", "CC_Mid_Anterior", "CC_Anterior"
)

# FS7 renamed Thalamus-Proper -> Thalamus; accept either
ASEG_ALIASES <- list(
  "Left.Thalamus.Proper"  = c("Left-Thalamus-Proper", "Left-Thalamus"),
  "Right.Thalamus.Proper" = c("Right-Thalamus-Proper", "Right-Thalamus")
)

FS_VERSIONS <- c("FS7_T1", "FS6_T1", "FS7_T1T2", "FS6_T1T2", "FS5.3",
                 "FS6_T1FLAIR", "FS7_T1FLAIR", "reconall-clinical", "FSInfant",
                 "synthseg", "synthseg-young", "synthseg-robust", "Custom")

dict_columns <- function() {
  regional <- unlist(lapply(c("lh", "rh"), function(h) {
    unlist(lapply(DK_REGIONS, function(r) paste(h, c("SA", "GM", "CT"), r, sep = ".")))
  }))
  c(regional,
    paste0("SUBC.", SUBC_STRUCTS),
    "fs_version_SA", "fs_version_CT", "fs_version_GM", "study_site",
    "GMV", "sGMV", "WMV", "CSF", "CBV", "TBV", "mean.CT", "total.SA",
    "sexMale", "logAge_days")
}

CHARACTER_COLUMNS <- c("fs_version_SA", "fs_version_CT", "fs_version_GM", "study_site")

# ---- argument parsing ------------------------------------------------------

DEFAULTS <- list(
  subjects_dir = NULL, subjects_list = NULL,
  lh_area = NULL, rh_area = NULL, lh_volume = NULL, rh_volume = NULL,
  lh_thickness = NULL, rh_thickness = NULL, aseg = NULL,
  demographics = NULL, demo_id_col = "participant",
  age_col = NULL, age_units = "days", age_type = "birth",
  ga_col = NULL, ga_units = "days",
  sex_col = NULL,
  male_values = "M,m,Male,male,MALE",
  female_values = "F,f,Female,female,FEMALE",
  study = NULL, site_col = NULL,
  fs_version = NULL, fs_version_sa = NULL, fs_version_ct = NULL,
  fs_version_gm = NULL,
  id_name = "participant",
  out = NULL
)

# printed by --help; the mapping details stay in the comment block at the top
USAGE <- c(
  "Reshape FreeSurfer output into the wide, one-row-per-scan CSV described by",
  "data_dictionary.csv (68 DK regions x {SA, GM, CT}, 30 subcortical volumes,",
  "global tissue volumes, pipeline version, site, sex, log age).",
  "",
  "Two input modes:",
  "",
  "1. recon-all subject directories -- reads",
  "     <subj>/stats/{lh.aparc.stats, rh.aparc.stats, aseg.stats}",
  "",
  "     Rscript scoring_new_data/freesurfer_to_datadict.R \\",
  "       --subjects-dir /path/to/SUBJECTS_DIR \\",
  "       --demographics demo.csv --demo-id-col participant \\",
  "       --age-col age_days --sex-col sex \\",
  "       --study ABCD --site-col site \\",
  "       --fs-version FS7_T1 \\",
  "       --out abcd_datadict.csv",
  "",
  "2. group tables from aparcstats2table / asegstats2table",
  "",
  "     Rscript scoring_new_data/freesurfer_to_datadict.R \\",
  "       --lh-area lh.area.tsv --rh-area rh.area.tsv \\",
  "       --lh-volume lh.vol.tsv --rh-volume rh.vol.tsv \\",
  "       --lh-thickness lh.thick.tsv --rh-thickness rh.thick.tsv \\",
  "       --aseg aseg.tsv \\",
  "       --demographics demo.csv --demo-id-col participant \\",
  "       --age-col age_years --age-units years --sex-col sex \\",
  "       --fs-version FS6_T1 \\",
  "       --out abcd_datadict.csv"
)

parse_args <- function(argv) {
  if (length(argv) == 0 || any(argv %in% c("-h", "--help"))) {
    cat(USAGE, sep = "\n")
    cat("\noptions:\n  --",
        paste(gsub("_", "-", names(DEFAULTS)), collapse = "\n  --"), "\n", sep = "")
    # don't take an interactive session down with us
    if (interactive()) stop("no arguments given; see the usage above", call. = FALSE)
    quit(status = 0)
  }
  opts <- DEFAULTS
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
    has_eq <- grepl("=", a, fixed = TRUE)
    flag <- sub("^--", "", if (has_eq) sub("=.*$", "", a) else a)
    key <- gsub("-", "_", flag)   # --subjects-dir and --subjects_dir both work
    if (!key %in% names(DEFAULTS)) stop("unknown option: --", flag, call. = FALSE)
    if (has_eq) {
      val <- sub("^[^=]*=", "", a)
    } else {
      if (i == length(argv)) stop("missing value for --", flag, call. = FALSE)
      val <- argv[i + 1]; i <- i + 1
    }
    opts[[key]] <- val
    i <- i + 1
  }
  opts
}

# ---- stats file parsers ----------------------------------------------------

read_stats_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hdr <- grep("^# ColHeaders", lines, value = TRUE)
  if (length(hdr) != 1L) stop("no unique '# ColHeaders' line in ", path, call. = FALSE)
  cols <- strsplit(trimws(sub("^# ColHeaders", "", hdr)), "[[:space:]]+")[[1]]
  body <- lines[!startsWith(lines, "#") & nzchar(trimws(lines))]
  if (!length(body)) stop("no data rows in ", path, call. = FALSE)
  utils::read.table(text = body, col.names = cols, stringsAsFactors = FALSE)
}

# "# Measure Cortex, CortexVol, Total cortical gray matter volume, 526511.2, mm^3"
# -> both "Cortex" and "CortexVol" key the value 526511.2
read_stats_measures <- function(path) {
  lines <- grep("^# Measure ", readLines(path, warn = FALSE), value = TRUE)
  out <- c()
  for (ln in lines) {
    f <- trimws(strsplit(sub("^# Measure[[:space:]]+", "", ln), ",")[[1]])
    num <- suppressWarnings(as.numeric(f))
    if (!any(!is.na(num))) next
    val <- num[max(which(!is.na(num)))]
    keys <- unique(f[seq_len(min(2L, length(f)))])
    out[keys] <- val
  }
  out
}

# named vector of canonical measurements for one subject's stats/ directory
read_subject_stats <- function(subj_dir) {
  vals <- c()
  for (hemi in c("lh", "rh")) {
    f <- file.path(subj_dir, "stats", paste0(hemi, ".aparc.stats"))
    if (!file.exists(f)) {
      warning("missing ", f, call. = FALSE)
      next
    }
    tab <- read_stats_table(f)
    aparc_cols <- c(SA = "SurfArea", GM = "GrayVol", CT = "ThickAvg")
    for (key in names(aparc_cols)) {
      meas <- aparc_cols[[key]]
      if (!meas %in% names(tab)) {
        warning("column ", meas, " not in ", f, call. = FALSE)
        next
      }
      v <- stats::setNames(tab[[meas]], paste(hemi, key, tab$StructName, sep = "."))
      vals[names(v)] <- v
    }
  }
  f <- file.path(subj_dir, "stats", "aseg.stats")
  if (file.exists(f)) {
    tab <- read_stats_table(f)
    v <- stats::setNames(tab$Volume_mm3, tab$StructName)
    vals[names(v)] <- v
    m <- read_stats_measures(f)
    vals[names(m)] <- m
  } else {
    warning("missing ", f, call. = FALSE)
  }
  vals
}

# ---- *stats2table readers --------------------------------------------------

read_group_table <- function(path) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(tab[[1]])
  mat <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- ids
  mat
}

# aparcstats2table columns look like "lh_bankssts_area" / "rh_insula_thickness"
aparc_table_vals <- function(path, hemi, key) {
  mat <- read_group_table(path)
  cn <- colnames(mat)
  region <- sub(paste0("^", hemi, "[_.]"), "", cn)
  region <- sub("_(area|volume|thickness|thicknessstd|meancurv)$", "", region)
  keep <- region %in% DK_REGIONS
  mat <- mat[, keep, drop = FALSE]
  colnames(mat) <- paste(hemi, key, region[keep], sep = ".")
  mat
}

# asegstats2table columns are StructNames plus global Measure columns
aseg_table_vals <- function(path) read_group_table(path)

# ---- assembling one row of data -----------------------------------------
# rename FS outputs to match data dictionary & add covariate placeholders.
# finds FS variable name that maps to dictionary var - pick first candidate name present in vals, else NA
pick <- function(vals, candidates) {
  hit <- candidates[candidates %in% names(vals)]
  if (!length(hit)) return(NA_real_)
  as.numeric(vals[[hit[1]]])
}

aseg_candidates <- function(struct) {
  alias <- ASEG_ALIASES[[struct]]
  unique(c(struct, gsub(".", "-", struct, fixed = TRUE), alias))
}

assemble_row <- function(vals) {
  row <- stats::setNames(rep(NA_real_, length(dict_columns())), dict_columns())

  #add hemi, SA/GM/CT string  
  for (hemi in c("lh", "rh")) {
    for (key in c("SA", "GM", "CT")) {
      for (r in DK_REGIONS) {
        nm <- paste(hemi, key, r, sep = ".")
        row[[nm]] <- pick(vals, nm)
      }
    }
  }

  for (s in SUBC_STRUCTS) {
    row[[paste0("SUBC.", s)]] <- pick(vals, aseg_candidates(s))
  }
  
  #map global phenos
  gmv  <- pick(vals, c("CortexVol", "Cortex"))
  sgmv <- pick(vals, c("SubCortGrayVol", "SubCortGray"))
  wmv  <- pick(vals, c("CerebralWhiteMatterVol", "CorticalWhiteMatterVol",
                       "CerebralWhiteMatter", "CorticalWhiteMatter"))
  bsv  <- pick(vals, c("BrainSegVol", "BrainSeg"))
  bsnv <- pick(vals, c("BrainSegVolNotVent", "BrainSegNotVent"))
  cbv  <- sum(row[paste0("SUBC.", c("Left.Cerebellum.Cortex",
                                    "Right.Cerebellum.Cortex",
                                    "Left.Cerebellum.White.Matter",
                                    "Right.Cerebellum.White.Matter"))])
  #rename/calc as needed and store in row
  row[["GMV"]]  <- gmv
  row[["sGMV"]] <- sgmv
  row[["WMV"]]  <- wmv
  row[["CSF"]]  <- bsv - bsnv
  row[["CBV"]]  <- cbv
  row[["TBV"]]  <- sum(c(wmv, sgmv, gmv, cbv))

  ct_cols <- grep("\\.CT\\.", dict_columns(), value = TRUE)
  sa_cols <- grep("\\.SA\\.", dict_columns(), value = TRUE)
  row[["mean.CT"]]  <- mean(row[ct_cols]) 
  row[["total.SA"]] <- sum(row[sa_cols])

  row
}

# ---- demographics ----------------------------------------------------------

split_csv_arg <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])

recode_sex <- function(x, male_values, female_values) {
  x <- trimws(as.character(x))
  out <- rep(NA_real_, length(x))
  out[x %in% male_values]   <- 1
  out[x %in% female_values] <- 0
  unmapped <- unique(x[is.na(out) & nzchar(x) & !is.na(x)])
  if (length(unmapped)) {
    warning("unmapped sex value(s) set to NA: ", paste(unmapped, collapse = ", "),
            call. = FALSE)
  }
  out
}

to_days <- function(x, units) {
  x <- as.numeric(x)
  switch(units,
         days = x,
         weeks = x * 7,
         months = x * (365.25 / 12),
         years = x * 365.25,
         stop("unsupported age units: ", units, call. = FALSE))
}

compute_log_age <- function(demo, opts) {
  if (is.null(opts$age_col)) {
    warning("no --age-col given; logAge_days left as NA", call. = FALSE)
    return(rep(NA_real_, nrow(demo)))
  }
  if (!opts$age_col %in% names(demo)) {
    stop("--age-col '", opts$age_col, "' not in demographics", call. = FALSE)
  }
  age_days <- to_days(demo[[opts$age_col]], opts$age_units)
  if (opts$age_type == "birth") {
    if (!is.null(opts$ga_col)) {
      if (!opts$ga_col %in% names(demo)) {
        stop("--ga-col '", opts$ga_col, "' not in demographics", call. = FALSE)
      }
      ga <- to_days(demo[[opts$ga_col]], opts$ga_units)
      ga[is.na(ga)] <- 280
    } else {
      ga <- rep(280, nrow(demo))
    }
    age_days <- age_days + ga
  } else if (opts$age_type != "postconception") {
    stop("--age-type must be 'birth' or 'postconception'", call. = FALSE)
  }
  if (any(age_days <= 0, na.rm = TRUE)) {
    warning("non-positive post-conception age(s); logAge_days set to NA there",
            call. = FALSE)
    age_days[age_days <= 0] <- NA_real_
  }
  log10(age_days)
}

# ---- main ------------------------------------------------------------------

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_args(argv)
  if (is.null(opts$out)) stop("--out is required", call. = FALSE)

  table_mode <- any(!vapply(opts[c("lh_area", "rh_area", "lh_volume", "rh_volume",
                                   "lh_thickness", "rh_thickness", "aseg")],
                            is.null, logical(1)))
  if (is.null(opts$subjects_dir) && !table_mode) {
    stop("give either --subjects-dir or the *stats2table files", call. = FALSE)
  }
  if (!is.null(opts$subjects_dir) && table_mode) {
    stop("use --subjects-dir or the *stats2table files, not both", call. = FALSE)
  }

  # ---- gather per-subject measurements
  if (!is.null(opts$subjects_dir)) {
    #just subjects listed
    if (!is.null(opts$subjects_list)) {
      subjects <- trimws(readLines(opts$subjects_list, warn = FALSE))
      subjects <- subjects[nzchar(subjects)]
    } else {
      #all subjects in dir
      subjects <- basename(dirname(dirname(
        Sys.glob(file.path(opts$subjects_dir, "*", "stats", "aseg.stats")))))
    }
    if (!length(subjects)) stop("no subjects with stats/aseg.stats under ",
                                opts$subjects_dir, call. = FALSE)
    message("found ", length(subjects), " subject(s)")
    rows <- lapply(subjects, function(s) {
      assemble_row(read_subject_stats(file.path(opts$subjects_dir, s)))
    })
    
    #table mode
  } else {
    tabs <- list()
    add <- function(m) {
      if (!length(tabs)) return(m)
      ids <- union(rownames(tabs), rownames(m))
      pad <- function(x) {
        out <- matrix(NA_real_, length(ids), ncol(x),
                      dimnames = list(ids, colnames(x)))
        out[rownames(x), ] <- x
        out
      }
      cbind(pad(tabs), pad(m))
    }
    specs <- list(c("lh_area", "lh", "SA"), c("rh_area", "rh", "SA"),
                  c("lh_volume", "lh", "GM"), c("rh_volume", "rh", "GM"),
                  c("lh_thickness", "lh", "CT"), c("rh_thickness", "rh", "CT"))
    for (sp in specs) {
      if (is.null(opts[[sp[1]]])) next
      tabs <- add(aparc_table_vals(opts[[sp[1]]], sp[2], sp[3]))
    }
    if (!is.null(opts$aseg)) tabs <- add(aseg_table_vals(opts$aseg))
    subjects <- rownames(tabs)
    message("found ", length(subjects), " subject(s) in the input tables")
    rows <- lapply(subjects, function(s) {
      v <- tabs[s, ]
      assemble_row(v[!is.na(v)])
    })
  }
  
  #stitch into one dataframe
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  df <- cbind(stats::setNames(data.frame(subjects, stringsAsFactors = FALSE),
                              opts$id_name), df)

  # ---- pipeline version
  fsv <- function(specific) if (!is.null(specific)) specific else opts$fs_version
  versions <- list(fs_version_SA = fsv(opts$fs_version_sa),
                   fs_version_CT = fsv(opts$fs_version_ct),
                   fs_version_GM = fsv(opts$fs_version_gm))
  for (v in names(versions)) {
    if (is.null(versions[[v]])) {
      warning(v, " not set (pass --fs-version); left as NA", call. = FALSE)
      df[[v]] <- NA_character_
    } else {
      if (!versions[[v]] %in% FS_VERSIONS) {
        warning(v, " value '", versions[[v]], "' is not one of the dictionary's ",
                "allowed levels: ", paste(FS_VERSIONS, collapse = ", "), call. = FALSE)
      }
      df[[v]] <- versions[[v]]
    }
  }

  # ---- demographics: study_site, sexMale, logAge_days
  if (!is.null(opts$demographics)) {
    demo <- utils::read.csv(opts$demographics, stringsAsFactors = FALSE,
                            colClasses = c(stats::setNames("character", opts$demo_id_col)))
    if (!opts$demo_id_col %in% names(demo)) {
      stop("--demo-id-col '", opts$demo_id_col, "' not in ", opts$demographics, call. = FALSE)
    }
    idx <- match(df[[opts$id_name]], trimws(demo[[opts$demo_id_col]]))
    unmatched <- df[[opts$id_name]][is.na(idx)]
    if (length(unmatched)) {
      warning(length(unmatched), " subject(s) had no demographics row (e.g. ",
              paste(utils::head(unmatched, 3), collapse = ", "), ")", call. = FALSE)
    }
    demo <- demo[idx, , drop = FALSE]

    if (!is.null(opts$sex_col)) {
      if (!opts$sex_col %in% names(demo)) {
        stop("--sex-col '", opts$sex_col, "' not in demographics", call. = FALSE)
      }
      df$sexMale <- recode_sex(demo[[opts$sex_col]],
                               split_csv_arg(opts$male_values),
                               split_csv_arg(opts$female_values))
    } else {
      warning("no --sex-col given; sexMale left as NA", call. = FALSE)
    }

    df$logAge_days <- compute_log_age(demo, opts)

    if (!is.null(opts$site_col)) {
      if (!opts$site_col %in% names(demo)) {
        stop("--site-col '", opts$site_col, "' not in demographics", call. = FALSE)
      }
      site <- as.character(demo[[opts$site_col]])
      site_col <- if (!is.null(opts$study)) paste(opts$study, site, sep = "_") else site
      site_col[is.na(site)] <- NA_character_   # paste() would give "STUDY_NA"
      df$study_site <- site_col
    } else {
      df$study_site <- NA_character_
    }
  } else {
    warning("no --demographics given; sexMale, logAge_days, and study_site left as NA",
            call. = FALSE)
  }

  # ---- final shape: dictionary order, dictionary classes
  keep <- c(opts$id_name, dict_columns())
  df <- df[, keep, drop = FALSE]
  for (v in intersect(CHARACTER_COLUMNS, names(df))) df[[v]] <- as.character(df[[v]])

  # ---- missingness, reported once for the whole run
  if (nrow(df) > 0) {
    vars <- setdiff(names(df), opts$id_name)
    n_na <- vapply(df[vars], function(x) sum(is.na(x)), integer(1))

    partial <- n_na[n_na > 0 & n_na < nrow(df)]
    if (length(partial)) {
      message("columns missing for some subjects (of ", nrow(df), " row(s)):")
      shown <- utils::head(sort(partial, decreasing = TRUE), 20)
      for (v in names(shown)) message("  ", v, ": ", shown[[v]], " NA")
      if (length(partial) > length(shown)) {
        message("  ... and ", length(partial) - length(shown), " more")
      }
    }

    # drop empty cols
    empty <- names(n_na)[n_na == nrow(df)]
    if (length(empty)) {
      message("dropping ", length(empty), " all-NA column(s): ",
              paste(empty, collapse = ", "))
      warning("dropped ", length(empty), " all-NA column(s): ",
              paste(empty, collapse = ", "), call. = FALSE)
      df <- df[, setdiff(names(df), empty), drop = FALSE]
    }
  }

  dir.create(dirname(opts$out), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, opts$out, row.names = FALSE, na = "")
  message("wrote ", nrow(df), " row(s) x ", ncol(df), " column(s) to ", opts$out)
  invisible(df)
}

# run only when this file is the script Rscript was pointed at -- source()ing it
# from an interactive session or from another script just defines the functions
sourced_from_elsewhere <- any(vapply(sys.frames(), function(e) !is.null(e$ofile),
                                     logical(1)))
if (!interactive() && !sourced_from_elsewhere) main()
