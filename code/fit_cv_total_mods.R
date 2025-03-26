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
l.name <- as.character(args[3])
fs <- as.character(args[4])
total <- as.character(args[5])
save_path <- as.character(args[6])
weight_pts <- as.logical(args[7])
log_scale <- as.logical(args[8])

#drop extra variables
df <- df %>%
  dplyr::select(all_of(c(pheno, fs, total,
                         "logAge_days", "sexMale", "study_site", "sexMale_x_logAge", "age_days"))) %>%
  na.omit()

#inverse-weight by age w/in sex (written w help from gpt)
if (weight_pts == TRUE){
  n_bins <- 50
  df$age_bin <- cut(df$age_days, breaks = n_bins, include.lowest = TRUE)
  
  
  df <- df %>%
    mutate(sex = as.factor(sexMale)) %>%
    group_by(sex, age_bin) %>%
    mutate(bin_count = n()) %>%
    group_by(sex) %>%
    mutate(
      observed_prob = bin_count / sum(bin_count),
      uniform_prob = 1 / n_bins,
      raw_weight = uniform_prob / observed_prob,
      weight = raw_weight * (n() / sum(raw_weight))
    ) %>%
    ungroup() %>%
    select(-bin_count, -observed_prob, -uniform_prob, - raw_weight)
  w <- "weighted"
} else{
  w <- "unweighted"
}

#log-scale pheno if necessary
if (log_scale == TRUE){
  pheno_sym <- sym(pheno)
  
  df <- df %>%
    mutate(!!pheno_sym := ifelse(!!pheno_sym==0, 1, !!pheno_sym)) %>% #replace 0 with 1
    mutate(!!pheno_sym := log(!!pheno_sym, base=10)) #transform
}

#define lambdas to be tested - FROM CONFIG FILE
if (l.name == "NULL"){
  l <- NULL
} else {
  l <- as.numeric(l.name)
}

#sim data ONCE for centile fan plotting
print("simulate data for plotting")
sim_df <- sim_data(df, "logAge_days", factor_var="sexMale", special_term = "sexMale_x_logAge = sexMale * logAge_days")

#loop over nu terms #NEED TO ADD START FROM 
nu_list <- list(int = "1",
                site = "study_site",
                sex = "sexMale",
                age = "logAge_days",
                sexAge = "sexMale + logAge_days",
                siteAge = "study_site + logAge_days",
                siteSex = "study_site + sexMale",
                siteAgeSex = "study_site + logAge_days + sexMale"
                )

#loop over fs moments
moment_list <- c("none", "mu", "both", "all")

#initialize empty lists
mod_list <- c()
results_df <- data.frame()
summary_df <- data.frame()

#FIT MODEL
for (fs_include in moment_list){
  
  for (total_include in moment_list){

    for (nu in nu_list){
    nu_name <- names(nu_list)[nu_list==nu]
    
    print(paste("fitting model with lambda =", l, 
                ", fs in", fs_include, 
                ", total in", total_include, 
                "and nu = ", nu_name, w))
    
    #FIT BASIC MODEL
    model <- gamlss_3lambda_etiv(pheno,
                            lambda=l,
                            total_var=total, total_moment=total_include,
                            fs_ver=fs, fs_moment=fs_include, 
                            fam="GG", 
                            weight=weight_pts,
                            nu_form=nu,
                            start.from = "mod_list[[1]]") #use first model as starting point
  
    #if model isn't fit, skip to next loop
    if (is.null(model)) {
      message("model fitting failed")
      next
    } else {
      message("model fit")
    }
    
    m_name <- paste(paste0("fs", fs_include), nu_name, paste0("total", total_include), sep="_")
    mod_list[[m_name]] <- model
    
    saveRDS(model, file=paste0(save_path, "/model_objs/", pheno, "_", w, "_lambda", l.name, "_", m_name, "_mod.rds"))
  
   
    #COMPILE
      #BIC & AIC
      tmp_df <- data.frame(
        "AIC" = model$aic,
        "BIC" = model$sbc,
        "lambda" = l.name,
        "pheno" = pheno,
        "fs_moment" = fs_include,
        "total_moment" = total_include,
        "weight" = w,
        "nu" = nu_name
      )
      summary_df <- rbind(summary_df, tmp_df)
    }
  }
}
expected <- length(nu_list)*length(moment_list)*length(moment_list)
print(paste(length(mod_list), "of", expected, "models fit"))

#SAVE CSVs
print("saving csvs")
fwrite(summary_df, file=paste0(save_path, "/model_sums/", pheno, "_", w, "_lambda", l.name, "_summary.csv"))

print("finding lowest BIC")
best_bic <- summary_df %>%
  arrange(BIC) %>%
  slice_head(n=1)

best_bic$m_name <- paste(paste0("fs", best_bic$fs_moment), 
                         best_bic$nu, 
                         paste0("total", best_bic$total_moment), 
                         sep="_")

print(best_bic$m_name)

best_mod <- mod_list[[best_bic$m_name]]

#RENAME BEST MOD
file.rename(paste0(save_path, "/model_objs/", best_bic$pheno, "_", best_bic$weight, "_lambda", best_bic$lambda, "_", best_bic$m_name, "_mod.rds"),
            paste0(save_path, "/model_objs/", best_bic$pheno, "_", best_bic$weight, "_lambda", best_bic$lambda, "_", best_bic$m_name, "_BestMod.rds"))

print("compiling stats")

#CENTILE FAN PLOT
print("creating centile fan plots")
fan_plot <- make_centile_fan(gamlssModel=best_mod, df=df, x_var="logAge_days", color_var="sexMale",
                             get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                             sim_data_list = sim_df,
                             show_points=FALSE
                             ) +
  ggtitle(paste(pheno, "\nsmoothed w/ lamda=", best_bic$lambda, ",", best_bic$m_name, w)) +
  xlab("log Age Days")

ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_", w, "_lambda", best_bic$lambda, "_", best_bic$m_name, ".png"), fan_plot)

fan_plot <- make_centile_fan(gamlssModel=best_mod, df=df, x_var=total, color_var="sexMale",
                             get_peaks=FALSE, desiredCentiles=c(0.05, 0.25, 0.5, 0.75, 0.95),
                             show_points=FALSE) +
  ggtitle(paste(pheno, "\nsmoothed w/ lamda=", best_bic$lambda, ",", best_bic$m_name, w)) +
  xlab(total)

ggsave(file=paste0(save_path, "/centile_plots/", pheno, "_", w, "_lambda", best_bic$lambda, "_", best_bic$m_name, "total.png"), fan_plot)


#WORM PLOT
print("creating worm plot")
wp <- wp.taki(xvar=df$logAge_days, resid=resid(best_mod), n.inter=8) +
  ggtitle(paste(pheno, "\nsmoothed w/ lambda=", l.name))
ggsave(file=paste0(save_path, "/worm_plots/", pheno, "_", w, "_lambda", best_bic$lambda, "_", best_bic$m_name, ".png"), wp)

#centiles
results_df <- cent_cdf(best_mod, df, plot=FALSE, group="sexMale")
results_df$lambda <- best_bic$lambda
results_df$fs <- best_bic$fs_moment
results_df$total <- best_bic$total_moment
results_df$nu <- best_bic$nu
resuts_df$weight <- w

#centiles
fwrite(results_df, file=paste0(save_path, "/cent_csvs/", pheno, "_", w, "_lambda", best_bic$lambda, "_", best_bic$m_name, "_results.csv"))

print("DONE")