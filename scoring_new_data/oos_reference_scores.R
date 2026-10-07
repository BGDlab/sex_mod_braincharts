# Script for obtaining reference (centile/z scores) from new data 
# benchmarked on brain charts from Gardner et al.

##################################
### LOAD PACKAGES
##################################

library(dplyr)
library(data.table)
library(gamlss)
#(re)install from github unless the installed copy came from `ref`.
#`ref` can be a branch name or a commit SHA
install_github_ref <- function(repo, ref) {
  pkg <- basename(repo)
  installed_ref <- if (requireNamespace(pkg, quietly = TRUE)) {
    utils::packageDescription(pkg)$RemoteRef
  }
  if (is.null(installed_ref) || !identical(installed_ref, ref)) {
    message("installing ", repo, "@", ref, " (found: ",
            if (is.null(installed_ref)) "none" else installed_ref, ")")
    remotes::install_github(paste0(repo, "@", ref), upgrade = "never",
                            build_vignettes = FALSE, force = TRUE)
  }
}

#--skip_install: packages were already checked (e.g. once by submit_oos_reference_scores.sh)
#--install_only: check/install packages, then exit
argv <- commandArgs(trailingOnly = TRUE)
skip_install <- "--skip_install" %in% argv
install_only <- "--install_only" %in% argv
argv <- argv[!argv %in% c("--skip_install", "--install_only")]

if (!skip_install) {
  install_github_ref("BGDlab/gamlssTools", "dev")
  install_github_ref("andy1764/gamlss2charts", "dev")
}
if (install_only) quit(save = "no", status = 0)
library(gamlssTools)
library(gamlss2charts)

##################################
### LOAD HELPER FUNCTIONS
##################################

#source helpers from github; if that fails (e.g. no internet), fall back to
#a local copy next to this script
helper_url <- "https://raw.githubusercontent.com/BGDlab/sex_mod_braincharts/refs/heads/sharing/scoring_new_data/oos_reference_scores_helper_funs.R"

script_dir <- function() {
  #Rscript
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) > 0) return(dirname(normalizePath(gsub("~+~", " ", f[1], fixed = TRUE))))
  #source()'d interactively
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  getwd()
}

tryCatch(
  source(helper_url),
  error = function(e) {
    local_helpers <- file.path(script_dir(), "oos_reference_scores_helper_funs.R")
    if (!file.exists(local_helpers))
      stop("could not source helper functions from GitHub (", conditionMessage(e),
           ") and no local copy found at ", local_helpers, call. = FALSE)
    warning("could not source helper functions from GitHub (", conditionMessage(e),
            "); using local copy at ", local_helpers, call. = FALSE)
    source(local_helpers)
  }
)

##################################
### DEFINE ARGUMENTS
##################################

#--pheno <name>: score this one phenotype instead of --pheno_list 
# for use with slurm job arrays
pheno <- NULL
i <- which(argv == "--pheno" | startsWith(argv, "--pheno="))
if (length(i) > 1) stop("--pheno given more than once", call. = FALSE)
if (length(i) == 1) {
  if (startsWith(argv[i], "--pheno=")) {
    pheno <- sub("^--pheno=", "", argv[i]); drop <- i
  } else {
    if (i == length(argv) || startsWith(argv[i + 1], "--"))
      stop("missing value for --pheno", call. = FALSE)
    pheno <- argv[i + 1]; drop <- c(i, i + 1)
  }
  if (!nzchar(pheno)) stop("empty value for --pheno", call. = FALSE)
  if (any(grepl("^--pheno[_-]list", argv)))
    stop("use --pheno or --pheno_list, not both", call. = FALSE)
  argv <- argv[-drop]
}

args <- parse_args(argv)
print(args)

#data to score
df <- fread(args$df, na.strings = c("NA", "", '""'))

# phenotypes to score - defaults to all
pheno_list <- if (!is.null(pheno)) pheno else readLines(args$pheno_list, warn = FALSE)

##################################
### VALIDATE INPUT DATA
##################################

#check variables used in all phenotype models
check_df(df)

df <- df %>%
  mutate(sexMale_x_logAge = sexMale * logAge_days, #calculate sex x age interaction
         .row_id = seq_len(dplyr::n())) #stable row key

#with a local model directory, fail now if it is missing or incomplete
check_model_dir(args$model_dir, pheno_list, args$total)

#when streaming, fail now if the models aren't on github at model_ref
check_model_remote(args$model_dir, pheno_list, args$total, args$model_ref)

##################################
### CALCULATE REFERENCE SCORES
##################################
print("calculating centiles...")

#score each pheno, holding its warnings and errors to print at the end of the run
results <- lapply(pheno_list, score_pheno_logged,
                  df = df,
                  total = args$total,
                  batch = args$batch,
                  ref_data = args$ref_data,
                  min_ref = args$min_ref,
                  model_dir = args$model_dir,
                  model_ref = args$model_ref)

df_cent <- Filter(Negate(is.null), lapply(results, `[[`, "scores"))
log_df  <- do.call(rbind, lapply(results, `[[`, "log"))

#rejoin to the full input data by row key
df_full_cent <- Reduce(
  function(a, b) dplyr::left_join(a, b, by = ".row_id"),
  df_cent,
  init = as.data.frame(df)
) %>%
  dplyr::select(-.row_id)

print(paste0("scored ", length(df_cent), "/", length(pheno_list), " phenotypes on ",
             nrow(df_full_cent), " subjects"))

##################################
### SAVE OUTPUTS
##################################
fwrite(df_full_cent, args$out_file)

##################################
### REPORT WARNINGS & ERRORS
##################################
#identical messages (e.g. package warnings repeated for every pheno) are shown once,
#with the phenotypes they came from
if (is.null(log_df)) {
  cat("\nno warnings or errors\n")
} else {
  for (type in c("error", "warning")) {
    l <- log_df[log_df$type == type, , drop = FALSE]
    if (nrow(l) == 0) next
    cat("\n==== ", toupper(type), "S (", nrow(l), " from ", length(unique(l$pheno)),
        " phenotype(s)) ====\n", sep = "")
    for (msg in unique(l$message)) {
      phenos <- unique(l$pheno[l$message == msg])
      cat("\n", msg, "\n", sep = "")
      if (length(phenos) > 1 || !grepl(phenos, msg, fixed = TRUE))
        cat("  [", length(phenos), " phenotype(s): ", paste(utils::head(phenos, 5), collapse = ", "),
            if (length(phenos) > 5) ", ...", "]\n", sep = "")
    }
  }
}
