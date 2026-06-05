#replot centiles, hardcoded for now
set.seed(99999)

#LOAD PACKAGES
library(data.table)
library(ggplot2)
library(dplyr)
library(gamlss)
library(gamlssTools)
library(cowplot)

base_path <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base_path, "code/gamlss_fit_funs.R"))
csv_path <- paste0(base_path, "cv_sample_?_test/*totalTRUE*")

#get pheno list
sa_list     <- readRDS(file.path(base_path, "pheno_lists/cortical_surf.rds"))
ct_list     <- readRDS(file.path(base_path, "pheno_lists/cortical_thickness.rds"))
vol_list    <- readRDS(file.path(base_path, "pheno_lists/cortical_vols.rds"))
global_list <- readRDS(file.path(base_path, "pheno_lists/global_vols.rds"))
sub_list    <- readRDS(file.path(base_path, "pheno_lists/subcortical_vols.rds"))
pheno_list  <- c(global_list, sub_list, vol_list, sa_list, ct_list)

pdf_fname <- paste0(base_path, "figs/trajectory_plts/all_phenos.pdf")
pdf(pdf_fname, width = 8, height = 10, bg = "white")

for (pheno in pheno_list) {
  f_list <- Sys.glob(paste0(csv_path, "/cent_csvs/", pheno, "_sexdiffs.csv"))
  f_list <- f_list[!grepl("weighted", f_list)] #ignore weighted models
  if (length(f_list) == 0) {
    warning(paste(length(f_list), "no files found for", pheno, "- skipping"))
    next
  } else if (length(f_list) == 1) {
    warning(paste(length(f_list), "1 file found for", pheno))
  }
  #get split
  names(f_list) <- sub(".*cv_sample_(.).*", "\\1", f_list)
  df <- rbindlist(lapply(f_list, fread), idcol = "split", use.names=TRUE)
  
 med_diff_plt <- ggplot(df) +
   geom_hline(yintercept=0, color="gray20", linetype = "dashed") +
    geom_line(aes(x=logAge_days, y=centile_M_minus_F_z, color=split)) +
   format_x_axis("log_lifespan_fetal", df$logAge_days) +
   labs(color= "Split-Half",
        y="Sex Difference in Median Trajectory",
        x= "Age at Scan (years)") +
   theme_minimal() +
   theme(legend.position = "none")
 
 var_diff_plt <- ggplot(df) +
   geom_hline(yintercept=0, color="gray20", linetype = "dashed") +
   geom_line(aes(x=logAge_days, y=logcv_M_div_F, color=split)) +
   format_x_axis("log_lifespan_fetal", df$logAge_days) +
   labs(color= "Split-Half",
        y= "Sex Difference in Variability Trajectory",
        x= "Age at Scan (years)") +
   theme_minimal() +
   theme(legend.position="bottom")
   
 
 shared_legend <- get_legend(var_diff_plt, "bottom")
 var_diff_plt <- var_diff_plt + theme(legend.position = "none")
 
 p <- plot_grid(
   med_diff_plt,
   var_diff_plt,
   shared_legend,
   ncol = 1,
   nrow = 3,
   rel_heights = c(1,1,.05)
 )
 
 #add a title so each pdf page is identifiable by pheno
 p <- plot_grid(
   ggdraw() + draw_label(pheno, fontface = "bold"),
   p,
   ncol = 1,
   rel_heights = c(0.05, 1)
 )
 
 fname <- paste0(base_path, "figs/trajectory_plts/", pheno,".png")
 ggsave(fname, p, width=8, height=10, dpi=300, bg = "white")
 
 print(p)
  
}
dev.off()
