#replot centiles, hardcoded for now
set.seed(99999)

#LOAD PACKAGES
library(data.table)
library(grid)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(gamlss)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
model <- readRDS(args[2]) #model to fit
save_path <- as.character(args[3])
total <- as.logical(args[4])
pheno <- as.character(args[5])
traintest <- as.character(args[6])
log_age <- TRUE #hard coding

birth <- log(280, base=10)
    
if (log_age ==TRUE){
    xvar <- "logAge_days"
    st <- "sexMale_x_logAge = sexMale * logAge_days"
} else {
    xvar <- "age_days"
    st <- "sexMale_x_age = sexMale * age_days"
}

if (total == TRUE){
  pt_show <- FALSE
}
    
model$call$data <- "df"
model$call$family <- "BCCG"
    
print(df)
print(model)
  
# Plot: predicted (line) + original (points)
    
p <- make_centile_fan(model, df, 
                      x_var=xvar, 
                      color_var="sexMale", 
                      desiredCentiles = c(0.5, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                      special_term = st,
                      remove_cent_effect="study_site",
                      show_points = pt_show,
                      remove_point_effect = "study_site",
                      x_scale=un_log) +
      theme_linedraw() +
      #theme(plot.title = element_blank()) +
      xlab("Age days") +
      scale_color_discrete(name = "Sex", labels = c("Female", "Male")) +
      guides(fill=FALSE) +
      geom_vline(xintercept=280)


#name - pheno, TRAIN or TEST, _plot.png
fname <- print(save_path, pheno, "_", traintest, "_plot.png")

print("saving plot")
ggsave(filename = fname, plot = p, width = 8, height = 5)
