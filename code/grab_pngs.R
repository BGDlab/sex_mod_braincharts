library(grid)
library(gridExtra)
library(png)
library(tools)

# Define root directory
root_dirs <- Sys.glob("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/cv_sample_?_test/*/centile_plots/")
#root_dirs <- Sys.glob("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts")

args <- commandArgs(trailingOnly=TRUE)
fname <- as.character(args[1])
pattern <- if (length(args) >= 2) args[2] else ".*\\.png$"  # Optional pattern filter

# Get all PNGs from matching directories - use find command for faster recursive search
message("Discovering PNG files...")
png_files <- unlist(lapply(root_dirs, function(dir) {
  # Use system find command which is often faster than list.files for deep directories
  cmd <- paste0("find '", dir, "' -type f -name '*.png' 2>/dev/null")
  system(cmd, intern = TRUE)
}))

# Apply pattern filter if provided (e.g., only centile plots, worm plots, etc.)
if (pattern != ".*\\.png$") {
  png_files <- png_files[grepl(pattern, basename(png_files))]
}

message(sprintf("Found %d PNG files", length(png_files)))

if (length(png_files) == 0) {
  stop("No PNG files found matching criteria")
}

# Output PDF path
output_pdf <- paste0("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/", fname, ".pdf")

# Create PDF
pdf(output_pdf, width = 8, height = 10)

# Pre-allocate layout viewport outside loop for efficiency
layout_vp <- viewport(layout = grid.layout(2, 1, heights = unit(c(1, 9), "null")))

# Process with progress tracking
total <- length(png_files)
success_count <- 0L
failed_files <- character(0)

for (i in seq_along(png_files)) {
  file <- png_files[i]
  
  # Progress message every 10% or every 50 files, whichever is more frequent
  if (i %% max(1, min(50, floor(total/10))) == 0 || i == 1) {
    message(sprintf("Processing %d/%d (%.1f%%)...", i, total, 100*i/total))
  }
  
  # Read the image
  img <- tryCatch({
    readPNG(file, native = FALSE)  # native=FALSE is faster for grid.raster
  }, error = function(e) {
    message(sprintf("Error reading %s: %s", basename(file), e$message))
    return(NULL)
  })
  
  if (!is.null(img)) {
    # Extract filename only (no path)
    filename <- basename(file)
    
    # Set up the layout (reuse pre-allocated viewport)
    grid.newpage()
    pushViewport(layout_vp)
    
    # Print filename
    grid.text(label = filename, 
              vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
              gp = gpar(fontsize = 14, fontface = "bold"))
    
    # Print image
    grid.raster(img, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    
    success_count <- success_count + 1L
  } else {
    failed_files <- c(failed_files, file)
  }
  
  # Clear image from memory immediately
  rm(img)
  gc(verbose = FALSE)  # Silent garbage collection
}

# Close PDF
dev.off()

message(sprintf("PDF saved as: %s", output_pdf))
message(sprintf("Successfully processed %d/%d images", success_count, total))
if (length(failed_files) > 0) {
  message(sprintf("Failed to process %d files", length(failed_files)))
}

