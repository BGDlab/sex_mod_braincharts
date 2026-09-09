#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# freesurfer_to_datadict.R
#
# Reshape FreeSurfer output into the wide, one-row-per-scan CSV described by
# data_dictionary.csv:
#   68 Desikan-Killiany regions x {SA, GM, CT}  (aparc.stats)
#   30 subcortical volumes  (aseg.stats)
#   global tissue volumes GMV / sGMV / WMV / CSF / CBV / TBV, mean.CT, total.SA
#   fs_version_SA / _CT / _GM, study_site, sexMale, logAge_days
#
# Two input modes:
#
# 1. recon-all subject directories (default) -- reads
#      <subj>/stats/{lh.aparc.stats, rh.aparc.stats, aseg.stats}
#
#      Rscript code/freesurfer_to_datadict.R \
#        --subjects-dir /path/to/SUBJECTS_DIR \
#        --demographics demo.csv --id-col participant \
#        --age-col age_days --sex-col sex \
#        --study ABCD --site-col site \
#        --fs-version FS7_T1 \
#        --out abcd_datadict.csv
#
# 2. group tables from aparcstats2table / asegstats2table
#
#      Rscript code/freesurfer_to_datadict.R \
#        --lh-area lh.area.tsv --rh-area rh.area.tsv \
#        --lh-volume lh.vol.tsv --rh-volume rh.vol.tsv \
#        --lh-thickness lh.thick.tsv --rh-thickness rh.thick.tsv \
#        --aseg aseg.tsv \
#        --demographics demo.csv --id-col participant \
#        --age-col age_years --age-units years --sex-col sex \
#        --study-site ABCD_site01 --fs-version FS6_T1 \
#        --out abcd_datadict.csv
#
# Demographics are joined on --id-col, matched against the FreeSurfer subject
# ID. Everything else is computed from the stats. Base R only, no packages.
#
# Mappings and assumptions:
#   SA = aparc SurfArea, GM = aparc GrayVol, CT = aparc ThickAvg
#   GMV  = CortexVol, sGMV = SubCortGrayVol,
#   WMV  = CerebralWhiteMatterVol (CorticalWhiteMatterVol in FS 5.3)
#   CSF  = BrainSegVol - BrainSegVolNotVent
#   CBV  = L/R cerebellum cortex + L/R cerebellum white matter
#   TBV  = WMV + sGMV + GMV + CBV
#   mean.CT / total.SA = unweighted mean / sum over the 68 DK regions; NA if any
#     region is missing, so partial parcellations are flagged rather than hidden
#   SUBC.*.Thalamus.Proper accepts the FS7 "Left-Thalamus" / "Right-Thalamus" name
#   logAge_days = log10(age_days + 280) unless --age-type postconception
#     (280 is overridable with --gestational-days / a per-subject --ga-col)
#   sexMale = 1/0 via --male-values / --female-values; anything else is NA
#   The data dictionary defines no ID variable, but the subject ID is written as
#     the first column so the file can be joined downstream; --no-id drops it.
#
# --dict data_dictionary.csv checks the output's names, order, and classes
# against the dictionary and reports any column containing NAs.
#
# Run with --help for the full option list.
# ---------------------------------------------------------------------------

# ---- column definitions (order matches data_dictionary.csv) ----------------

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
# string with "." -> "-", except where ASEG_ALIASES says otherwise
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
  `subjects-dir` = NULL, `subjects-list` = NULL,
  `lh-area` = NULL, `rh-area` = NULL, `lh-volume` = NULL, `rh-volume` = NULL,
  `lh-thickness` = NULL, `rh-thickness` = NULL, aseg = NULL,
  demographics = NULL, `id-col` = NULL,
  `age-col` = NULL, `age-units` = "days", `age-type` = "birth",
  `ga-col` = NULL, `ga-units` = "days", `gestational-days` = "280",
  `sex-col` = NULL,
  `male-values` = "1,M,m,Male,male,MALE",
  `female-values` = "0,F,f,Female,female,FEMALE",
  study = NULL, `site-col` = NULL, `study-site` = NULL,
  `fs-version` = NULL, `fs-version-sa` = NULL, `fs-version-ct` = NULL,
  `fs-version-gm` = NULL,
  `id-name` = "participant", `no-id` = FALSE,
  dict = NULL, out = NULL
)

parse_args <- function(argv) {
  if (length(argv) == 0 || any(argv %in% c("-h", "--help"))) {
    header <- readLines(script_path(), warn = FALSE)
    header <- header[seq_len(grep("^# ---- column definitions", header)[1] - 1)]
    cat(sub("^#", "", header[startsWith(header, "#") & !startsWith(header, "#!")]),
        sep = "\n")
    cat("\noptions:\n  --", paste(names(DEFAULTS), collapse = "\n  --"), "\n", sep = "")
    quit(status = 0)
  }
  opts <- DEFAULTS
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
    if (grepl("=", a, fixed = TRUE)) {
      key <- sub("^--", "", sub("=.*$", "", a))
      val <- sub("^[^=]*=", "", a)
    } else {
      key <- sub("^--", "", a)
      if (key == "no-id") { val <- TRUE } else {
        if (i == length(argv)) stop("missing value for --", key, call. = FALSE)
        val <- argv[i + 1]; i <- i + 1
      }
    }
    if (!key %in% names(DEFAULTS)) stop("unknown option: --", key, call. = FALSE)
    opts[[key]] <- val
    i <- i + 1
  }
  opts
}

script_path <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) normalizePath(f) else "code/freesurfer_to_datadict.R"
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

# ---- assembling one dictionary row -----------------------------------------

pick <- function(vals, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(vals)]
  if (!length(hit)) {
    if (required) warning("could not find ", label, " (looked for: ",
                          paste(candidates, collapse = ", "), ")", call. = FALSE)
    return(NA_real_)
  }
  as.numeric(vals[[hit[1]]])
}

aseg_candidates <- function(struct) {
  alias <- ASEG_ALIASES[[struct]]
  unique(c(struct, gsub(".", "-", struct, fixed = TRUE), alias))
}

assemble_row <- function(vals, id) {
  row <- stats::setNames(rep(NA_real_, length(dict_columns())), dict_columns())

  for (hemi in c("lh", "rh")) {
    for (key in c("SA", "GM", "CT")) {
      for (r in DK_REGIONS) {
        nm <- paste(hemi, key, r, sep = ".")
        row[[nm]] <- pick(vals, nm, paste(id, nm))
      }
    }
  }

  for (s in SUBC_STRUCTS) {
    row[[paste0("SUBC.", s)]] <- pick(vals, aseg_candidates(s), paste(id, s))
  }

  gmv  <- pick(vals, c("CortexVol", "Cortex", "TotalGrayVol"), paste(id, "GMV"))
  sgmv <- pick(vals, c("SubCortGrayVol", "SubCortGray"), paste(id, "sGMV"))
  wmv  <- pick(vals, c("CerebralWhiteMatterVol", "CorticalWhiteMatterVol",
                       "CerebralWhiteMatter", "CorticalWhiteMatter"),
               paste(id, "WMV"))
  bsv  <- pick(vals, c("BrainSegVol", "BrainSeg"), paste(id, "BrainSegVol"))
  bsnv <- pick(vals, c("BrainSegVolNotVent", "BrainSegNotVent"),
               paste(id, "BrainSegVolNotVent"))
  cbv  <- sum(row[paste0("SUBC.", c("Left.Cerebellum.Cortex",
                                    "Right.Cerebellum.Cortex",
                                    "Left.Cerebellum.White.Matter",
                                    "Right.Cerebellum.White.Matter"))])

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
         months = x * (365.245 / 12),
         years = x * 365.245,
         stop("unsupported age units: ", units, call. = FALSE))
}

compute_log_age <- function(demo, opts) {
  if (is.null(opts$`age-col`)) {
    warning("no --age-col given; logAge_days left as NA", call. = FALSE)
    return(rep(NA_real_, nrow(demo)))
  }
  if (!opts$`age-col` %in% names(demo)) {
    stop("--age-col '", opts$`age-col`, "' not in demographics", call. = FALSE)
  }
  age_days <- to_days(demo[[opts$`age-col`]], opts$`age-units`)
  if (opts$`age-type` == "birth") {
    if (!is.null(opts$`ga-col`)) {
      if (!opts$`ga-col` %in% names(demo)) {
        stop("--ga-col '", opts$`ga-col`, "' not in demographics", call. = FALSE)
      }
      ga <- to_days(demo[[opts$`ga-col`]], opts$`ga-units`)
      ga[is.na(ga)] <- as.numeric(opts$`gestational-days`)
    } else {
      ga <- rep(as.numeric(opts$`gestational-days`), nrow(demo))
    }
    age_days <- age_days + ga
  } else if (opts$`age-type` != "postconception") {
    stop("--age-type must be 'birth' or 'postconception'", call. = FALSE)
  }
  if (any(age_days <= 0, na.rm = TRUE)) {
    warning("non-positive post-conception age(s); logAge_days set to NA there",
            call. = FALSE)
    age_days[age_days <= 0] <- NA_real_
  }
  log10(age_days)
}

# ---- validation ------------------------------------------------------------

validate_output <- function(df, opts) {
  expected <- dict_columns()
  if (!is.null(opts$dict)) {
    dict <- utils::read.csv(opts$dict, stringsAsFactors = FALSE)
    expected <- dict$variable
    got <- setdiff(names(df), c(if (!isTRUE(opts$`no-id`)) opts$`id-name`))
    if (!identical(got, expected)) {
      missing <- setdiff(expected, got)
      extra   <- setdiff(got, expected)
      if (length(missing)) message("!! missing dictionary variables: ",
                                   paste(missing, collapse = ", "))
      if (length(extra)) message("!! variables not in dictionary: ",
                                 paste(extra, collapse = ", "))
      if (!length(missing) && !length(extra)) message("!! column order differs from dictionary")
    } else {
      message("dictionary check: all ", length(expected),
              " variables present, in dictionary order")
    }
    bad_class <- character(0)
    for (v in intersect(expected, names(df))) {
      want <- dict$class[match(v, dict$variable)]
      is_num <- is.numeric(df[[v]])
      if (want == "numeric" && !is_num) bad_class <- c(bad_class, v)
      if (want == "character" && is_num) bad_class <- c(bad_class, v)
    }
    if (length(bad_class)) message("!! class mismatch: ",
                                   paste(bad_class, collapse = ", "))
  }
  na_counts <- vapply(df[intersect(expected, names(df))],
                      function(x) sum(is.na(x)), integer(1))
  incomplete <- na_counts[na_counts > 0]
  if (length(incomplete)) {
    message("columns with missing values (n rows = ", nrow(df), "):")
    for (v in names(incomplete)) message("  ", v, ": ", incomplete[[v]], " NA")
  } else {
    message("no missing values")
  }
  invisible(NULL)
}

# ---- main ------------------------------------------------------------------

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_args(argv)
  if (is.null(opts$out)) stop("--out is required", call. = FALSE)

  table_mode <- any(!vapply(opts[c("lh-area", "rh-area", "lh-volume", "rh-volume",
                                   "lh-thickness", "rh-thickness", "aseg")],
                            is.null, logical(1)))
  if (is.null(opts$`subjects-dir`) && !table_mode) {
    stop("give either --subjects-dir or the *stats2table files", call. = FALSE)
  }
  if (!is.null(opts$`subjects-dir`) && table_mode) {
    stop("use --subjects-dir or the *stats2table files, not both", call. = FALSE)
  }

  # ---- gather per-subject measurements
  if (!is.null(opts$`subjects-dir`)) {
    if (!is.null(opts$`subjects-list`)) {
      subjects <- trimws(readLines(opts$`subjects-list`, warn = FALSE))
      subjects <- subjects[nzchar(subjects)]
    } else {
      subjects <- basename(dirname(dirname(
        Sys.glob(file.path(opts$`subjects-dir`, "*", "stats", "aseg.stats")))))
    }
    if (!length(subjects)) stop("no subjects with stats/aseg.stats under ",
                                opts$`subjects-dir`, call. = FALSE)
    message("found ", length(subjects), " subject(s)")
    rows <- lapply(subjects, function(s) {
      assemble_row(read_subject_stats(file.path(opts$`subjects-dir`, s)), s)
    })
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
    specs <- list(c("lh-area", "lh", "SA"), c("rh-area", "rh", "SA"),
                  c("lh-volume", "lh", "GM"), c("rh-volume", "rh", "GM"),
                  c("lh-thickness", "lh", "CT"), c("rh-thickness", "rh", "CT"))
    for (sp in specs) {
      if (is.null(opts[[sp[1]]])) next
      tabs <- add(aparc_table_vals(opts[[sp[1]]], sp[2], sp[3]))
    }
    if (!is.null(opts$aseg)) tabs <- add(aseg_table_vals(opts$aseg))
    subjects <- rownames(tabs)
    message("found ", length(subjects), " subject(s) in the input tables")
    rows <- lapply(subjects, function(s) {
      v <- tabs[s, ]
      assemble_row(v[!is.na(v)], s)
    })
  }

  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  df <- cbind(stats::setNames(data.frame(subjects, stringsAsFactors = FALSE),
                              opts$`id-name`), df)

  # ---- pipeline version
  fsv <- function(specific) if (!is.null(specific)) specific else opts$`fs-version`
  versions <- list(fs_version_SA = fsv(opts$`fs-version-sa`),
                   fs_version_CT = fsv(opts$`fs-version-ct`),
                   fs_version_GM = fsv(opts$`fs-version-gm`))
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
    if (is.null(opts$`id-col`)) stop("--id-col is required with --demographics",
                                     call. = FALSE)
    demo <- utils::read.csv(opts$demographics, stringsAsFactors = FALSE,
                            colClasses = c(stats::setNames("character", opts$`id-col`)))
    if (!opts$`id-col` %in% names(demo)) {
      stop("--id-col '", opts$`id-col`, "' not in ", opts$demographics, call. = FALSE)
    }
    idx <- match(df[[opts$`id-name`]], trimws(demo[[opts$`id-col`]]))
    unmatched <- df[[opts$`id-name`]][is.na(idx)]
    if (length(unmatched)) {
      warning(length(unmatched), " subject(s) had no demographics row (e.g. ",
              paste(utils::head(unmatched, 3), collapse = ", "), ")", call. = FALSE)
    }
    demo <- demo[idx, , drop = FALSE]

    if (!is.null(opts$`sex-col`)) {
      if (!opts$`sex-col` %in% names(demo)) {
        stop("--sex-col '", opts$`sex-col`, "' not in demographics", call. = FALSE)
      }
      df$sexMale <- recode_sex(demo[[opts$`sex-col`]],
                               split_csv_arg(opts$`male-values`),
                               split_csv_arg(opts$`female-values`))
    } else {
      warning("no --sex-col given; sexMale left as NA", call. = FALSE)
    }

    df$logAge_days <- compute_log_age(demo, opts)

    if (!is.null(opts$`site-col`)) {
      if (!opts$`site-col` %in% names(demo)) {
        stop("--site-col '", opts$`site-col`, "' not in demographics", call. = FALSE)
      }
      site <- as.character(demo[[opts$`site-col`]])
      df$study_site <- if (!is.null(opts$study)) paste(opts$study, site, sep = "_") else site
    }
  } else {
    warning("no --demographics given; sexMale and logAge_days left as NA",
            call. = FALSE)
  }

  if (is.null(opts$`site-col`)) {
    if (!is.null(opts$`study-site`)) {
      df$study_site <- opts$`study-site`
    } else if (!is.null(opts$study)) {
      df$study_site <- opts$study
    } else {
      warning("no --study-site / --study / --site-col; study_site left as NA",
              call. = FALSE)
      df$study_site <- NA_character_
    }
  }

  # ---- final shape: dictionary order, dictionary classes
  keep <- c(if (!isTRUE(opts$`no-id`)) opts$`id-name`, dict_columns())
  df <- df[, keep, drop = FALSE]
  for (v in intersect(CHARACTER_COLUMNS, names(df))) df[[v]] <- as.character(df[[v]])

  validate_output(df, opts)

  dir.create(dirname(opts$out), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, opts$out, row.names = FALSE, na = "")
  message("wrote ", nrow(df), " row(s) x ", ncol(df), " column(s) to ", opts$out)
  invisible(df)
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  if (!interactive()) main()
}
