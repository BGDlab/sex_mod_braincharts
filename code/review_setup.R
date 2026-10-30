### shared setup for review_test_models.Rmd & review_trajectories.Rmd ###
# split out of review_test_models.Rmd (setup chunk + integral funs)

if(!require('pacman')) {
  install.packages('pacman')
}
pacman::p_load(data.table, ggplot2, tidyverse, stats, skimr, rlang, ggrepel, gamlss, ggseg, ggsegExtra, RESI, paletteer, patchwork, plotly, pracma, gcplyr, legendry, gghighlight, forcats, sf, cowplot, subcortexVisualizationR, rstatix, freesurferformats, mgcv)

#devtools::install_github("BGDlab/gamlssTools", force=TRUE)
library(gamlssTools)

sa_list <- readRDS("./pheno_lists/cortical_surf.rds")
ct_list <- readRDS("./pheno_lists/cortical_thickness.rds")
vol_list <- readRDS("./pheno_lists/cortical_vols.rds")
global_list <- readRDS("./pheno_lists/global_vols.rds")
sub_list <- readRDS("./pheno_lists/subcortical_vols.rds")
pheno_list <- c(sa_list, ct_list, vol_list, global_list, sub_list)
test_path <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_revision/cv_sample_?_test*/*/"
#visualizing ventricles separately
vent_list <- c("CSF", "SUBC.Left.Lateral.Ventricle", "SUBC.Right.Lateral.Ventricle", "SUBC.3rd.Ventricle", "SUBC.4th.Ventricle")

# install.packages("remotes")
#remotes::install_github("BGDlab/ggsegTissue", force=TRUE)
library(ggsegTissue)

#source spin test
source("https://raw.githubusercontent.com/frantisekvasa/rotate_parcellation/refs/heads/master/R/perm.sphere.p.R")
source("https://raw.githubusercontent.com/frantisekvasa/rotate_parcellation/refs/heads/master/R/rotate.parcellation.R")
fs_path <- "/mnt/isilon/bgdlab_processing/code/freesurfer-6.0/subjects/fsaverage/"

source("code/plot_cv_brain.R")

### integral funs ###

#fun for derivative
get_deriv <- function(x, age) {
  c(NA, diff(x) / diff(age))
}

#riemann_stieltjes - updating default fun based on claude/convo with taki
rs_integral <- function(x, y, g = function(x) 10^(x)) {
  stopifnot(
    length(x) == length(y),
    length(x) >= 2
  )

  # sort just in case
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  gx <- g(x)

  # monotonicity check (important!)
  if (any(diff(gx) < 0)) {
    stop("g(x) must be non-decreasing for a Riemann–Stieltjes integral")
  }

  sum(y[-length(y)] * diff(gx))
}

#euclidean distance fun from https://www.geeksforgeeks.org/r-language/how-to-calculate-euclidean-distance-in-r/
euclidean_dist <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))


#filter out rows indicating sensitivity analysis, except where passed to keep
sensitivity_cols <- c("weighted", "age2plus", "runnerup")

filter_sensitivity <- function(df, keep = NULL, sens = sensitivity_cols) {
  df %>%
    dplyr::filter(dplyr::if_all(dplyr::any_of(setdiff(sens, keep)), ~ !as.logical(.x)),
                  dplyr::if_all(dplyr::any_of(keep), ~ as.logical(.x)))
}