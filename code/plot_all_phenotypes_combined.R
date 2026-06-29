# Combine train-vs-test phenotype plots and sex-difference trajectory plots
# into a single PDF: all train/test pages first, then all sex-diff trajectory pages.
set.seed(99999)

library(grid)
library(png)
library(data.table)
library(ggplot2)
library(dplyr)
library(gamlss)
library(gamlssTools)
library(cowplot)

base_dir <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts"
source(file.path(base_dir, "code/gamlss_fit_funs.R"))

# get pheno list
sa_list     <- readRDS(file.path(base_dir, "pheno_lists/cortical_surf.rds"))
ct_list     <- readRDS(file.path(base_dir, "pheno_lists/cortical_thickness.rds"))
vol_list    <- readRDS(file.path(base_dir, "pheno_lists/cortical_vols.rds"))
global_list <- readRDS(file.path(base_dir, "pheno_lists/global_vols.rds"))
sub_list    <- readRDS(file.path(base_dir, "pheno_lists/subcortical_vols.rds"))
pheno_list  <- c(global_list, sub_list, vol_list, sa_list, ct_list)

output_pdf <- file.path(base_dir, "all_phenotypes_combined.pdf")
pdf(output_pdf, width = 8, height = 10, bg = "white")

## ---------------------------------------------------------------------------
## Part 1: train-vs-test phenotype plots
## ---------------------------------------------------------------------------

directions <- list(
  list(train = "A", test = "B"),
  list(train = "B", test = "A")
)
total_levels <- c("TRUE", "FALSE")

# Glob across category subdirs (global_vols/cortical_*/subcortical_vols) to find
# the appropriate <pheno>_plot.png for a given fold sample, role, and total flag.
find_png <- function(sample, role, total, pheno) {
  pattern <- file.path(
    base_dir,
    sprintf("cv_sample_%s_%s", sample, role),
    sprintf("*_total%s*logAgeTRUE_pbmods", total),
    "replot",
    sprintf("%s_*_plot.png", pheno)
  )
  hits <- Sys.glob(pattern)
  hits <- hits[!grepl("weighted", hits, fixed = TRUE)]
  if (length(hits) == 0) NA_character_ else hits[1]
}

# Header band over two equal-height image panels
layout_vp <- viewport(layout = grid.layout(3, 1, heights = unit(c(1, 9, 9), "null")))

total_pages <- length(pheno_list) * length(directions) * length(total_levels)
page_idx       <- 0L
written        <- 0L
missing_pages  <- character(0)

for (pheno in pheno_list) {
  for (dir in directions) {
    for (tot in total_levels) {
      page_idx <- page_idx + 1L
      if (page_idx == 1L || page_idx %% 50L == 0L) {
        message(sprintf("Page %d/%d (%s, train=%s/test=%s, total=%s)",
                        page_idx, total_pages, pheno, dir$train, dir$test, tot))
      }

      top_path    <- find_png(dir$train, "train", tot, pheno)
      bottom_path <- find_png(dir$test,  "test",  tot, pheno)

      top_img <- if (!is.na(top_path)) {
        tryCatch(readPNG(top_path, native = FALSE), error = function(e) NULL)
      } else NULL
      bottom_img <- if (!is.na(bottom_path)) {
        tryCatch(readPNG(bottom_path, native = FALSE), error = function(e) NULL)
      } else NULL

      top_status    <- if (is.na(top_path))    "MISSING" else if (is.null(top_img))    "READ ERROR" else "ok"
      bottom_status <- if (is.na(bottom_path)) "MISSING" else if (is.null(bottom_img)) "READ ERROR" else "ok"

      if (top_status != "ok" || bottom_status != "ok") {
        missing_pages <- c(missing_pages, sprintf(
          "%s | train=%s test=%s total=%s | top=%s bottom=%s",
          pheno, dir$train, dir$test, tot, top_status, bottom_status
        ))
      }

      # Skip the page only when neither image is available.
      if (is.null(top_img) && is.null(bottom_img)) next

      grid.newpage()
      pushViewport(layout_vp)

      header <- sprintf("%s  |  train=%s / test=%s  |  total=%s",
                        pheno, dir$train, dir$test, tot)
      grid.text(header,
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
                gp = gpar(fontsize = 14, fontface = "bold"))

      if (!is.null(top_img)) {
        grid.raster(top_img,
                    vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
      } else {
        grid.text(sprintf("(train plot %s for cv_sample_%s_train)", top_status, dir$train),
                  vp = viewport(layout.pos.row = 2, layout.pos.col = 1),
                  gp = gpar(fontsize = 12, fontface = "italic", col = "grey40"))
      }
      if (!is.null(bottom_img)) {
        grid.raster(bottom_img,
                    vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      } else {
        grid.text(sprintf("(test plot %s for cv_sample_%s_test)", bottom_status, dir$test),
                  vp = viewport(layout.pos.row = 3, layout.pos.col = 1),
                  gp = gpar(fontsize = 12, fontface = "italic", col = "grey40"))
      }

      written <- written + 1L
      rm(top_img, bottom_img)
      gc(verbose = FALSE)
    }
  }
}

message(sprintf("Part 1 (train vs test): wrote %d of %d pages", written, total_pages))
if (length(missing_pages) > 0) {
  message(sprintf("Skipped %d pages due to missing/unreadable files:", length(missing_pages)))
  for (m in missing_pages) message("  ", m)
}

## ---------------------------------------------------------------------------
## Part 2: sex-difference trajectory plots
## ---------------------------------------------------------------------------

csv_path <- file.path(base_dir, "cv_sample_?_test/*totalTRUE*")

# y-limits spanning at least [-0.1, 0.1], widened to fit the data if needed.
min_ylim <- function(x) c(min(x, -0.1, na.rm = TRUE), max(x, 0.1, na.rm = TRUE))

for (pheno in pheno_list) {
  f_list <- Sys.glob(file.path(csv_path, "cent_csvs", sprintf("%s_sexdiffs.csv", pheno)))
  f_list <- f_list[!grepl("weighted", f_list)] # ignore weighted models
  if (length(f_list) == 0) {
    warning(paste(length(f_list), "no files found for", pheno, "- skipping"))
    next
  } else if (length(f_list) == 1) {
    warning(paste(length(f_list), "1 file found for", pheno))
  }
  # get split
  names(f_list) <- sub(".*cv_sample_(.).*", "\\1", f_list)
  df <- rbindlist(lapply(f_list, fread), idcol = "split", use.names = TRUE)

  med_diff_plt <- ggplot(df) +
    geom_hline(yintercept = 0, color = "gray20", linetype = "dashed") +
    geom_line(aes(x = logAge_days, y = centile_M_minus_F_z, color = split)) +
    format_x_axis("log_lifespan_fetal", df$logAge_days) +
    labs(color = "Split-Half",
         y = "Sex Difference in Median Trajectory",
         x = "Age at Scan (years)") +
    theme_minimal() +
    theme(legend.position = "none") +
    coord_cartesian(ylim = min_ylim(df$centile_M_minus_F_z))

  var_diff_plt <- ggplot(df) +
    geom_hline(yintercept = 0, color = "gray20", linetype = "dashed") +
    geom_line(aes(x = logAge_days, y = logcv_M_div_F, color = split)) +
    format_x_axis("log_lifespan_fetal", df$logAge_days) +
    labs(color = "Split-Half",
         y = "Sex Difference in Variability Trajectory",
         x = "Age at Scan (years)") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = min_ylim(df$logcv_M_div_F))

  shared_legend <- get_legend(var_diff_plt, "bottom")
  var_diff_plt <- var_diff_plt + theme(legend.position = "none")

  p <- plot_grid(
    med_diff_plt,
    var_diff_plt,
    shared_legend,
    ncol = 1,
    nrow = 3,
    rel_heights = c(1, 1, .05)
  )

  # add a title so each pdf page is identifiable by pheno
  p <- plot_grid(
    ggdraw() + draw_label(pheno, fontface = "bold"),
    p,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )

  fname <- file.path(base_dir, "figs/trajectory_plts", sprintf("%s.png", pheno))
  ggsave(fname, p, width = 8, height = 10, dpi = 300, bg = "white")

  print(p)
}

dev.off()

message(sprintf("PDF saved as: %s", output_pdf))
