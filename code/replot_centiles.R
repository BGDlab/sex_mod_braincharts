#replot centiles, hardcoded for now
set.seed(99999)

#LOAD PACKAGES
library(data.table)
library(grid)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(gamlss)
devtools::install_github("BGDlab/gamlssTools", force=TRUE, upgrade=FALSE)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))

# paths
# paths
dir_paths <- list()
dir_paths[[1]] <- list(
  csv_dir = "data/pheno_dfs_totalFALSE",
  rds_dir = "braincharts/cortical_surf_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods",
  csv_path = "totalFALSE_logPhenoTRUE_logAgeTRUE\\.csv$",
  log_age=TRUE
)
dir_paths[[2]] <- list(
  csv_dir = "data/pheno_dfs_totalFALSE",
  rds_dir = "braincharts/cortical_thickness_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods",
  csv_path = "totalFALSE_logPhenoTRUE_logAgeTRUE\\.csv$",
  log_age=TRUE
)
dir_paths[[3]] <- list(
  csv_dir = "data/pheno_dfs_totalFALSE",
  rds_dir = "braincharts/cortical_vols_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods",
  csv_path = "totalFALSE_logPhenoTRUE_logAgeTRUE\\.csv$",
  log_age=TRUE
)
dir_paths[[4]] <- list(
  csv_dir = "data/pheno_dfs_totalFALSE",
  rds_dir = "braincharts/global_vols_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods",
  csv_path = "totalFALSE_logPhenoTRUE_logAgeTRUE\\.csv$",
  log_age=TRUE
)
dir_paths[[5]] <- list(
  csv_dir = "data/pheno_dfs_totalFALSE",
  rds_dir = "braincharts/subcortical_vols_totalTRUE_logPhenoTRUE_logAgeTRUE_pbmods",
  csv_path = "totalTRUE_logPhenoTRUE_logAgeTRUE\\.csv$",
  log_age=TRUE
)


birth <- log(280, base=10)

unscale <- function(x){10^x - 5}

pdf(paste0(base,"braincharts_centiles_unscaledAge_july25.pdf"), width = 8.5, height = 11)  # Open PDF device

for (paths in dir_paths){
  # List all .rds files ending with BestMod.rds
  rds_dir <- paste0(paths$rds_dir,"/model_objs")
  csv_dir <- paths$csv_dir
  csv_path <- paths$csv_path
  log_age <- paths$log_age
  
  plot_dir <- paste0(paths$rds_dir,"/unscaled_centile_plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  rds_files <- list.files(rds_dir, pattern = "BestMod\\.rds$", full.names = TRUE)
  
  # List all .csv files
  csv_files <- list.files(csv_dir, pattern = csv_path, full.names = TRUE)
  
  # Extract matching prefixes from .rds filenames
  get_prefix <- function(path) sub("_.*$", "", basename(path))
  
  rds_prefixes <- setNames(rds_files, vapply(rds_files, get_prefix, ""))
  
  # Match each rds file with the corresponding csv file
  matched <- lapply(names(rds_prefixes), function(prefix) {
    # Find the CSV with the same prefix
    csv_path <- grep(paste0("^", prefix), basename(csv_files), value = TRUE)
    
    if (length(csv_path) == 1) {
      list(
        prefix = prefix,
        model = readRDS(rds_prefixes[[prefix]]),
        data  = fread(file.path(csv_dir, csv_path), stringsAsFactors=TRUE, na.strings = ""),
        log_age = log_age
      )
    } else {
      warning(sprintf("No unique match found for prefix: %s", prefix))
      NULL
    }
  })
  
  # Remove NULLs (if any unmatched)
  matched <- Filter(Negate(is.null), matched)
  
  for (item in matched) {
    model <- item$model
    df  <- item$data
    prefix <- item$prefix
    log_age <- item$log_age
    
    if (log_age ==TRUE){
      xvar <- "logAge_days"
      st <- "sexMale_x_logAge = sexMale * logAge_days"
    } else {
      xvar <- "age_days"
      st <- "sexMale_x_age = sexMale * age_days"
    }
    
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
      theme(plot.title = element_blank()) +
      xlab("Age (log days)") +
      scale_color_discrete(name = "Sex", labels = c("Female", "Male")) +
      guides(fill=FALSE) +
      geom_vline(xintercept=birth)
    
    model$call$data <- "df"
    model$call$family <- "BCCG"
    
    
    # results_df <- cent_cdf(model, df, plot=FALSE, group="sexMale")
    # 
    # tbl <- tableGrob(results_df)
    
    print(p)
    grid.newpage()
    # grid.draw(tbl)
    ggsave(filename = file.path(plot_dir, paste0(prefix, "_plot.png")), plot = p, width = 8, height = 5)
    
  }
}
dev.off()