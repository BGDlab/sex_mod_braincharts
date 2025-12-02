#replot centiles, hardcoded for now
set.seed(99999)

#LOAD PACKAGES
library(data.table)
library(grid)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(gamlss)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))

# paths
# paths
dir_paths <- list()
# dir_paths[[1]] <- list(
#   csv_dir = "data/pheno_dfs_totalFALSE",
#   rds_dir = "braincharts/cortical_surf_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods",
#   csv_path = "totalFALSE_logPhenoTRUE_logAgeTRUE\\.csv$",
#   log_age=TRUE
# )
# dir_paths[[2]] <- list(
#   csv_dir = "data/pheno_dfs_totalFALSE",
#   rds_dir = "braincharts/cortical_thickness_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods",
#   csv_path = "totalFALSE_logPhenoTRUE_logAgeTRUE\\.csv$",
#   log_age=TRUE
# )
# dir_paths[[3]] <- list(
#   csv_dir = "data/pheno_dfs_totalFALSE",
#   rds_dir = "braincharts/cortical_vols_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods",
#   csv_path = "totalFALSE_logPhenoTRUE_logAgeTRUE\\.csv$",
#   log_age=TRUE
# )
# dir_paths[[4]] <- list(
#   csv_dir = "data/pheno_dfs_totalFALSE",
#   rds_dir = "braincharts/global_vols_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods",
#   csv_path = "totalFALSE_logPhenoTRUE_logAgeTRUE\\.csv$",
#   log_age=TRUE
# )
dir_paths[[1]] <- list(
  csv_dir = "data/pheno_dfs_totalFALSE",
  rds_dir = "braincharts/subcortical_vols_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods",
  csv_path = "totalFALSE.*logAgeTRUE\\.csv$",  # Handle optional _logPheno*_ in filename
  log_age=TRUE
)

birth <- log(280, base=10)

unscale <- function(x){10^x - 5}

for (paths in dir_paths){
  # List all .rds files ending with BestMod.rds
  rds_dir <- paste0(paths$rds_dir,"/model_objs")
  csv_dir <- paths$csv_dir
  csv_path <- paths$csv_path
  log_age <- paths$log_age
  
  plot_dir <- paste0(paths$rds_dir,"/replots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE) #get best mod or full testing mod
  
  # List all .csv files
  csv_files <- list.files(csv_dir, pattern = csv_path, full.names = TRUE)
  
  # Extract matching prefixes from .rds filenames
  get_prefix <- function(path) sub("_.*$", "", basename(path))
  
  rds_prefixes <- setNames(rds_files, vapply(rds_files, get_prefix, ""))
  
  # Match each rds file with the corresponding csv file
  matched <- lapply(names(rds_prefixes), function(pr) {
    # Find the CSV with the same prefix
    csv_path2 <- grep(paste0("^", pr), basename(csv_files), value = TRUE)
    
    if (length(csv_path2) == 1) {
      list(
        prefix = pr,
        model = readRDS(rds_prefixes[[pr]]),
        dat = fread(file.path(csv_dir, csv_path2), stringsAsFactors=TRUE, na.strings = ""),
        log_age = log_age
      )
    } else {
      warning(sprintf("No unique match found for prefix: %s", pr))
      NULL
    }
  })
  
  # Remove NULLs (if any unmatched)
  matched <- Filter(Negate(is.null), matched)
  
  for (item in matched) {
    model <- item$model
    df  <- item$dat
    prefix <- item$prefix
    log_age <- item$log_age
    
    if (log_age ==TRUE){
      xvar <- "logAge_days"
      st <- "sexMale_x_logAge = sexMale * logAge_days"
    } else {
      xvar <- "age_days"
      st <- "sexMale_x_age = sexMale * age_days"
    }
    
    model$call$data <- "df"
    model$call$family <- "BCCG"
    
    print(df)
    print(model)
    # Plot: predicted (line) + original (points)
    p <- make_centile_fan(model, df, 
                          x_var=xvar, 
                          color_var="sexMale", 
                          desiredCentiles = c(0.5, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                          special_term = st,
                          remove_cent_effect="study_site",
                          remove_point_effect = "study_site",
                          y_scale=unscale,
                          x_scale=un_log) +
      theme_linedraw() +
      #theme(plot.title = element_blank()) +
      xlab("Age days") +
      scale_color_discrete(name = "Sex", labels = c("Female", "Male")) +
      guides(fill=FALSE) +
      geom_vline(xintercept=280)
    
    print(p)
    print("saving plot")
    ggsave(filename = file.path(plot_dir, paste0(prefix, "_plot.png")), plot = p, width = 8, height = 5)
    
  }
}
