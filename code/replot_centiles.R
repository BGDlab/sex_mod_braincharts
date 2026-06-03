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
total <- as.logical(args[4])
save_path <- as.character(args[5])

if (total == TRUE){
  pt_show <- FALSE
} else {
  pt_show <- TRUE
}
    
model$call$data <- "df"
model$call$family <- model$family[[1]]
    
print(df)
print(model)
  
#CENTILE FAN PLOT
pheno <- model$mu.terms[[2]] %>% as.character()

pred_list <- list_predictors(model)

# pick age var
if ("logAge_days" %in% pred_list){
  age_var      <- "logAge_days"
  sex_age_term <- "sexMale_x_logAge = sexMale * logAge_days"
  vars_of_interest <- c("sexMale_x_logAge", age_var, "sexMale")
  birth <- log(280, base=10)
} else {
  age_var      <- "age_days"
  sex_age_term <- "sexMale_x_age = sexMale * age_days"
  vars_of_interest <- c("sexMale_x_age", age_var, "sexMale")
  birth <- 280
}

# build special_term; if total-corrected, hold the total var at 0 for plotting
special_term <- sex_age_term
if (isTRUE(total)){
  total_candidates <- c("TBV", "mean.CT", "total.SA")
  total_var <- intersect(total_candidates, pred_list)
  if (length(total_var) == 1){
    special_term <- c(sex_age_term, paste0(total_var, "=0"))
    message("total==TRUE, zeroing ", total_var, " in sim_data")
  } else if (length(total_var) == 0){
    warning("total==TRUE but no known total var (TBV/mean.CT/total.SA) in pred_list; not zeroing")
  } else {
    stop("multiple total vars in pred_list: ", paste(total_var, collapse = ", "))
  }
}

print("simulate data for plotting")
sim_df <- sim_data(df, age_var, factor_var = "sexMale", special_term = special_term)

# residualize points
resid_terms <- setdiff(pred_list, vars_of_interest)
    
p <- make_centile_fan(gamlssModel=model, 
                      df=df, 
                      x_var=age_var, 
                      color_var="sexMale",
                      get_peaks=TRUE, 
                      desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                      sim_data_list = sim_df,
                      remove_point_effect = resid_terms,
                      x_axis="log_lifespan_fetal") +
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
ggsave(filename = fname, plot = p, width = 8, height = 5)
saveRDS(p, fname_rds)