#Fit GAMLSS models to select from on CV samples

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

source("./code/gamlss_fit_funs.R")

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
print(args)
df <- fread(args[1], stringsAsFactors = TRUE, na.strings = "") #path to csv
pheno <- as.character(args[2])
freesurf <- as.character(args[3])
save_path <- as.character(args[4])

#drop extra variables
df <- df %>%
  dplyr::select(all_of(c(pheno, freesurf, "logAge_days", "sexMale", "study_site", "sexMale_x_logAge", "age_days"))) %>%
  na.omit() %>%
  rename(freesurfer = freesurf) %>%
  trunc_coverage(vars=("logAge_days")) #drop points at ends if too sparse

#try just fitting MOST COMPLEX 4-param model and see what happens

mu_str <- "pb(logAge_days, control = pb.control(order = 3)) + pb(sexMale_x_logAge, control = pb.control(order = 3)) + sexMale + random(study_site) + freesurfer"
sig_str <- "pb(logAge_days, control = pb.control(order = 3)) + pb(sexMale_x_logAge, control = pb.control(order = 3)) + sexMale + random(study_site) + freesurfer"
nu_str <- "pb(logAge_days, control = pb.control(order = 3)) + pb(sexMale_x_logAge, control = pb.control(order = 3)) + sexMale + study_site + freesurfer"
tau_str <- "logAge_days + sexMale_x_logAge + sexMale + random(study_site) + freesurfer"
  
print("fitting model")
model <- gamlss_4param(pheno, 
                       mu_base=mu_str,
                       sig_base=sig_str,
                       nu_base=nu_str,
                       tau_base=tau_str,
                       fam="BCT",
                       start.from=NULL)

saveRDS(model, file=paste0(save_path, "/model_objs/", pheno, "_4param_mod.rds"))


#sim data ONCE for centile fan plotting
print("simulate data for plotting")
sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")

#CENTILE FAN PLOT
print("creating centile fan plot")
fan_plot <- make_centile_fan(gamlssModel=model, df=df, x_var="logAge_days", color_var="sexMale",
                             get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                             sim_data_list = sim_df)  +
  labs(title=paste(pheno, "validation model"),
       x ="log Age (days)",
       color = "Sex=Male", fill="Sex=Male")

ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_4param.png"), fan_plot)

#WORM PLOT
print("creating worm plot")
wp <- wp.taki(xvar=df$logAge_days, resid=resid(model), n.inter=8) +
  ggtitle(paste(pheno, "validation model"))
ggsave(file=paste0(save_path, "/worm_plots/", pheno, "_4param.png"), wp)

#COMPILE
print("compiling stats")
#centiles
results_df <- cent_cdf(model, df, "sexMale")

#BIC & AIC
summary_df <- data.frame(
  "AIC" = model$aic,
  "BIC" = model$sbc,
  "pheno" = pheno
)

#SAVE CSVs
print("saving csvs")

#centiles
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", pheno, "_4param_centiles.csv"))

#BIC & AIC
fwrite(summary_df, file=paste0(save_path, "/model_sums/", filename, "_4param_summary.csv"))
