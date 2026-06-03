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
traintest <- as.character(args[3])
split <- as.character(args[4])
total <- as.logical(args[5])
save_path <- as.character(args[6])

birth <- log(280, base=10)

if (total == TRUE){
  pt_show <- FALSE
} else {
  pt_show <- TRUE
}
    
model$call$data <- "df"
model$call$family <- "BCCG"
    
print(df)
print(model)
  
#CENTILE FAN PLOT
pheno <- model$mu.terms[[2]] %>% as.character()

pred_list <- list_predictors(model)

#check if age is log-scaled
if ("logAge_days" %in% pred_list){
  age_var <- "logAge_days"
  print("simulate data for plotting")
  sim_df <- sim_data(df, age_var, factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")
  vars_of_interest <- c("sexMale_x_logAge", age_var, "sexMale")
} else {
  print("simulate data for plotting")
  age_var <- "age_days"
  sim_df <- sim_data(df, age_var, factor_var="sexMale", special_term = "sexMale_x_age = sexMale * age_days")
  vars_of_interest <- c("sexMale_x_age", age_var, "sexMale")
}

#define nuisance covars to residualize from points
resid_terms <- setdiff(pred_list, vars_of_interest)
    
p <- make_centile_fan(model, df, 
                      x_var=age_var, 
                      color_var="sexMale", 
                      desiredCentiles = c(0.5, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                      sim_data_list = sim_df,
                      show_points = pt_show,
                      remove_point_effect = resid_terms,
                      x_scale=un_log) +
      theme_linedraw() +
      #theme(plot.title = element_blank()) +
      xlab("Age days") +
      scale_color_discrete(name = "Sex", labels = c("Female", "Male")) +
      guides(fill=FALSE) +
      geom_vline(xintercept=birth)


#name - pheno, TRAIN or TEST, _plot.png
fname     <- file.path(save_path, paste0(pheno, "_", traintest, "_plot.png"))
fname_rds <- file.path(save_path, paste0(pheno, "_", traintest, "_plot.rds"))

print(paste("saving plot to", fname))
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
ggsave(filename = fname, plot = p, width = 8, height = 5)
saveRDS(p, fname_rds)