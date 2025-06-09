library(grid)
library(gridExtra)
library(png)
library(tools)

# Define root directory
root_dirs <- Sys.glob("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/figs/") #cortical/subcort training mods
#root_dirs <- Sys.glob("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts")

# Get all PNGs from matching directories
png_files <- unlist(lapply(root_dirs, function(dir) {
  list.files(dir, pattern = "\\.png$", full.names = TRUE
             #, recursive = TRUE
             )
}))
print(length(png_files))

# Output PDF path
output_pdf <- "samp_w_slip.pdf"

# Create PDF
pdf(output_pdf, width = 8, height = 10)

for (file in png_files) {
  # Read the image
  img <- tryCatch({
    readPNG(file)
  }, error = function(e) {
    message(paste("Error reading", file))
    return(NULL)
  })
  
  if (!is.null(img)) {
    # Extract filename only (no path)
    filename <- basename(file)
    #print(filename)

    # Set up the layout
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 9), "null"))))
    
    # Print filename
    grid.text(label = filename, vp = viewport(layout.pos.row = 1, layout.pos.col = 1), gp = gpar(fontsize = 14, fontface = "bold"))
    
    # Print image
    grid.raster(img, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  }
}

# Close PDF
dev.off()

message("PDF saved as: ", output_pdf)

