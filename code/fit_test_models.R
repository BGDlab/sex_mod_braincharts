#sub-divide IDPs into the largest possible complete chunks for ComBat harmonization

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
pheno <- as.character(args[2])
save_path <- as.character(args[3])

#drop extra variables
df <- df %>%
  dplyr::select(all_of(c(pheno, "logAge_days", "sexMale", "fs_version", "study"))) %>%
  na.omit()

#define degrees of freedom to be tested
degree_list <- seq(3, 24, by=3)

results_df <- data.frame("degree" = as.numeric(),
                         "sigma_degree" = as.numeric(),
                         "pheno" = as.character(),
                         "BIC" = as.numeric(),
                         "AIC" = as.numeric()
                         )

#define gamlss fitting function
gamlss_try <- function(pheno, df, degree, sigma_degree){
  result <- tryCatch({
    eval(parse(text = paste0("gamlss(formula =", pheno, "~ ns(logAge_days, df = ", degree, ") + sexMale + fs_version + study,
                  sigma.formula = ~ ns(logAge_days, df = ", sigma_degree, ") + sexMale + fs_version + study,
                  nu.formula = ~ 1,
                  control = gamlss.control(n.cyc = 200), 
                  family = GG, data=", df, ", trace = FALSE)")))
  } , warning = function(w) {
    message("warning")
    eval(parse(text = paste0("gamlss(formula =", pheno, "~ ns(logAge_days, df = ", degree, ") + sexMale + fs_version + study,
                  sigma.formula = ~ ns(logAge_days, df = ", sigma_degree, ") + sexMale + fs_version + study,
                  nu.formula = ~ 1,
                  control = gamlss.control(n.cyc = 200), 
                  family = GG, data=", df, ", trace = FALSE)")))
  } , error = function(e) {
    message("error, trying method=CG()")
    eval(parse(text = paste0("gamlss(formula =", pheno, "~ ns(logAge_days, df = ", degree, ") + sexMale + fs_version + study,
                  sigma.formula = ~ ns(logAge_days, df = ", sigma_degree, ") + sexMale + fs_version + study,
                  nu.formula = ~ 1, method=CG(),
                  control = gamlss.control(n.cyc = 200), 
                  family = GG, data=", df, ", trace = FALSE)")))
  } , finally = {
    message("done")
  } )
}

#sim data ONCE for centile fan plotting
sim_df <- sim_data(df, "logAge_days", color_var= "sexMale")

#FIT BASE MODEL
for (degree in degree_list){
  sigma_degree <- degree/3
  
  model <- gamlss_try(pheno, df, degree, sigma_degree)
  
  saveRDS(model, file=paste0(save_path, "/model_objs/", pheno, "_", "mu", degree, "sig", sigma_degree, "mod.rds"))
  
  #save centile fan plot
  fan_plot <- make_centile_fan(gamlssModel=model, df=df, x_var="logAge_days", color_var="sexMale",
                               get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                               sim_data_list = sim_df) +
    ggtitle(paste(pheno, "smoothed\nw mu.df=", degree, ", sigma.df=", sigma_degree))
  
  ggsave(file=paste0(save_path, "/plots/", pheno, "_", "mu", degree, "sig", sigma_degree, ".png"), fan_plot)
  
  #compile results
  sub_df <- data.frame("degree" = degree,
                       "sigma_degree" = sigma_degree,
                       "pheno" = pheno,
                       "BIC" = BIC(model),
                       "AIC" = AIC(model))
  
  results_df <- rbind(results_df, sub_df)
}

fwrite(results_df, file=paste0(save_path, "/", pheno, "_results.csv"))

print()
