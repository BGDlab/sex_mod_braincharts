# Compile train-vs-test phenotype plots into one PDF.
# Each page = one (pheno x direction x totalTRUE/FALSE) combination,
# with the training-model plot on top and the testing-model plot below.
# Splits: direction A->B uses cv_sample_A_train (top) + cv_sample_B_test (bottom);
#         direction B->A uses cv_sample_B_train (top) + cv_sample_A_test (bottom).

library(grid)
library(png)

base_dir <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts"

sa_list     <- readRDS(file.path(base_dir, "pheno_lists/cortical_surf.rds"))
ct_list     <- readRDS(file.path(base_dir, "pheno_lists/cortical_thickness.rds"))
vol_list    <- readRDS(file.path(base_dir, "pheno_lists/cortical_vols.rds"))
global_list <- readRDS(file.path(base_dir, "pheno_lists/global_vols.rds"))
sub_list    <- readRDS(file.path(base_dir, "pheno_lists/subcortical_vols.rds"))
pheno_list  <- c(global_list, sub_list, vol_list, sa_list, ct_list)

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
  #print(pattern)
  hits <- Sys.glob(pattern)
  hits <- hits[!grepl("weighted", hits, fixed = TRUE)]
  if (length(hits) == 0) NA_character_ else hits[1]
}

output_pdf <- file.path(base_dir, "all_phenotypes_train_vs_test.pdf")
pdf(output_pdf, width = 8, height = 10)

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

dev.off()

message(sprintf("PDF saved as: %s", output_pdf))
message(sprintf("Wrote %d of %d pages", written, total_pages))
if (length(missing_pages) > 0) {
  message(sprintf("Skipped %d pages due to missing/unreadable files:", length(missing_pages)))
  for (m in missing_pages) message("  ", m)
}
