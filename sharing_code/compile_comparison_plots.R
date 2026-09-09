library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(png)
library(grid)

# usage: Rscript compile_comparison_plots.R <fig_path> [pheno_lists_dir] [out_filename]
#   fig_path        dir containing the *_zdiffs.rds files (and where the pdf is saved)
#   pheno_lists_dir dir of phenotype-category .rds lists (default "./pheno_lists")
#   out_filename    pdf filename, saved inside fig_path (default "zdiffs_histogram.pdf")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
fig_path        <- as.character(args[1])
pheno_lists_dir <- if (length(args) >= 2) as.character(args[2]) else "./pheno_lists"
out_filename    <- if (length(args) >= 3) as.character(args[3]) else "sharing_mod_comparison.pdf"

#gather z-diff rds files (all.files=TRUE so dot-prefixed/hidden files are included too)
zdiff_files <- list.files(fig_path, pattern = "_zdiffs\\.rds$", all.files = TRUE, full.names = TRUE)
if (length(zdiff_files) == 0) {
  stop("No *_zdiffs.rds files found in ", fig_path)
}
z_df <- map_dfr(zdiff_files, readRDS)

#build phenotype -> category lookup from pheno_lists dir
list_files <- list.files(pheno_lists_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(list_files) == 0) {
  stop("No phenotype-category .rds lists found in ", pheno_lists_dir)
}
pheno_lists <- setNames(lapply(list_files, readRDS),
                         tools::file_path_sans_ext(basename(list_files)))

pretty_labels <- c(
  global_vols         = "Global Vol",
  cortical_vols        = "Regional Vol",
  subcortical_vols     = "Subcortical Vol",
  cortical_surf        = "Regional SA",
  cortical_thickness   = "Regional CT"
)

get_pheno_cat <- function(pheno) {
  hits <- names(pheno_lists)[map_lgl(pheno_lists, ~ pheno %in% .x)]
  if (length(hits) == 0) NA_character_ else hits[1]
}

z_df <- z_df %>%
  mutate(pheno_cat = map_chr(phenotype, get_pheno_cat),
         pheno_cat = recode(pheno_cat, !!!pretty_labels))

#plot
plt_df <- z_df %>%
  pivot_longer(c(z_max, z_mean), names_to = "stat", values_to = "diff") %>%
  mutate(stat = recode(stat, z_max = "max_z", z_mean = "mean_z"))

plt <- ggplot(plt_df) +
  geom_histogram(aes(x = diff, fill = pheno_cat), position = "identity", alpha = .5) +
  facet_grid(total ~ stat) +
  scale_x_log10() +
  labs(x = "|difference| (log scale)") +
  theme_bw()

#gather already-saved centile fan pngs (all.files=TRUE so dot-prefixed/hidden files are included too)
centfan_files <- sort(list.files(fig_path, pattern = "_centfan\\.png$", all.files = TRUE, full.names = TRUE))

#write histogram + centile fans as pages of one pdf
#centile fan pngs are all saved at ggplot's default 7x7in (square), so a
#square page fits every image exactly with no per-image aspect check needed
pdf_path <- file.path(fig_path, out_filename)
pdf(pdf_path, width = 7, height = 7)
print(plt)
for (f in centfan_files) {
  img <- readPNG(f)
  grid.newpage()
  grid.raster(img)
}
dev.off()

message("saved histogram + ", length(centfan_files), " centile fan(s) to ", pdf_path)
