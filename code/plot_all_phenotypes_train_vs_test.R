# Compile all phenotype plots from train and test into one large plot using patchwork
# Each row is a phenotype: train (left), test (right)

library(patchwork)
library(ggplot2)

# Directories (edit as needed)
train_dir <- "cv_sample_A_train/global_vols_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods/replots"
test_dir  <- "cv_sample_B_test/global_vols_totalFALSE_logPhenoTRUE_logAgeTRUE_pbmods/replots"

# List all plot files
train_plots <- list.files(train_dir, pattern = "_plot\\.rds$", full.names = TRUE)
test_plots  <- list.files(test_dir,  pattern = "_plot\\.rds$", full.names = TRUE)

# Extract phenotype prefix from filename
get_prefix <- function(path) sub("_plot\\.rds$", "", basename(path))

train_prefixes <- setNames(train_plots, vapply(train_plots, get_prefix, ""))
test_prefixes  <- setNames(test_plots,  vapply(test_plots,  get_prefix, ""))

# Find common phenotypes
common_phenos <- intersect(names(train_prefixes), names(test_prefixes))

# For each phenotype, create a row (train | test)
plot_rows <- lapply(common_phenos, function(pheno) {
  train_plot <- readRDS(train_prefixes[[pheno]]) + ggtitle(paste(pheno, "Train"))
  test_plot  <- readRDS(test_prefixes[[pheno]])  + ggtitle(paste(pheno, "Test"))
  train_plot + test_plot
})

# Combine all rows into one big plot (2 columns, n rows)
big_plot <- wrap_plots(plot_rows, ncol = 1) +
  plot_annotation(title = "All Phenotypes: Train (left) vs Test (right)")

# Save
out_file <- "all_phenotypes_train_vs_test.png"
ggsave(out_file, plot = big_plot, width = 16, height = 5 * length(plot_rows))
cat("Saved combined plot to", out_file, "\n") 