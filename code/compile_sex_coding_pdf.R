#Compile refit_sex_coding.R outputs into one pdf: a page per BestMod comparing centile fans & sex-specific
#worm plots across sexMale (original), sexFemale, and sexEffect codings
#usage: Rscript compile_sex_coding_pdf.R <list of BestMod paths> <output.pdf>

library(grid)
library(gridExtra)
library(png)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
mod_list <- args[1] #same list passed to subjobs_sex_coding.sh: BestMod path [tab csv path] per line
output_pdf <- args[2]

mods <- readLines(mod_list)
mods <- trimws(sub("\t.*$", "", mods))
mods <- mods[mods != "" & !startsWith(mods, "#")]

codings <- c(female_ref = "sexMale\n(F=0, M=1)\noriginal",
             male_ref = "sexFemale\n(M=0, F=1)",
             effect = "sexEffect\n(F=-1, M=1)")

#draw a png into the current layout cell, or a note if it's missing
draw_png <- function(file, row, col) {
  vp <- viewport(layout.pos.row = row, layout.pos.col = col)
  if (file.exists(file)) {
    grid.raster(readPNG(file), vp = vp)
  } else {
    grid.text("not available", vp = vp, gp = gpar(col = "grey50"))
  }
}

pdf(output_pdf, width = 17, height = 13)

n_missing <- 0
for (i in seq_along(mods)) {
  mod_file <- mods[i]
  m_name <- sub("_BestMod\\.rds$", "", basename(mod_file))
  out_dir <- file.path(dirname(dirname(mod_file)), "sex_coding", m_name)
  message(sprintf("%d/%d: %s", i, length(mods), out_dir))

  stats_file <- file.path(out_dir, "sex_coding_stats.csv")
  stats <- if (file.exists(stats_file)) fread(stats_file) else NULL
  if (is.null(stats)) n_missing <- n_missing + 1

  title <- paste0(basename(dirname(dirname(dirname(mod_file)))), " / ",
                  basename(dirname(dirname(mod_file))), "\n", m_name)

  grid.newpage()
  #rows: title, header, 3 codings, stats table; cols: coding label, fan, female worm, male worm
  pushViewport(viewport(layout = grid.layout(
    6, 4,
    heights = unit(c(0.7, 0.3, 3, 3, 3, 1.5), "null"),
    widths = unit(c(0.6, 4, 2.4, 2.4), "null"))))
  grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4),
            gp = gpar(fontsize = 16, fontface = "bold"))
  headers <- c("Centile fan", "Female worm plot", "Male worm plot")
  for (j in 1:3) {
    grid.text(headers[j], vp = viewport(layout.pos.row = 2, layout.pos.col = j + 1),
              gp = gpar(fontsize = 13, fontface = "bold"))
  }

  for (k in seq_along(codings)) {
    cd <- names(codings)[k]
    grid.text(codings[[cd]], vp = viewport(layout.pos.row = k + 2, layout.pos.col = 1),
              gp = gpar(fontsize = 12, fontface = "bold"))
    draw_png(file.path(out_dir, paste0("fan_", cd, ".png")), k + 2, 2)
    draw_png(file.path(out_dir, paste0("worm_", cd, "_F.png")), k + 2, 3)
    draw_png(file.path(out_dir, paste0("worm_", cd, "_M.png")), k + 2, 4)
  }

  #fit stats across codings
  if (!is.null(stats)) {
    keep <- intersect(c("sex_var", "converged", "deviance", "df_fit", "AIC", "BIC", "mu_sexAge_edf", "sigma_sexAge_edf"), names(stats))
    tab <- as.data.frame(stats)[, keep, drop = FALSE]
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], function(x) formatC(x, format = "f", digits = 2))
    tab$n <- stats$n
    pushViewport(viewport(layout.pos.row = 6, layout.pos.col = 1:4))
    grid.draw(tableGrob(tab, rows = NULL, theme = ttheme_minimal(base_size = 11)))
    popViewport()
  } else {
    grid.text(paste("no outputs found in", out_dir), vp = viewport(layout.pos.row = 6, layout.pos.col = 1:4),
              gp = gpar(col = "firebrick"))
  }
  popViewport()
}

dev.off()
message("PDF saved as: ", output_pdf)
if (n_missing > 0) message(n_missing, " of ", length(mods), " models had no refit outputs")
