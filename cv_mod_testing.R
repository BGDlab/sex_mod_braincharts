library(dplyr)
library(gamlss)
library(gamlssTools)
library(data.table)

pheno <- "WMV"

df <- fread("./data/cv_sample_A.csv", stringsAsFactors = TRUE)

df <- df %>%
  dplyr::select(all_of(c(pheno, "fs_version_GM", "logAge_days", "sexMale", "study_site", "sexMale_x_logAge", "age_days"))) %>%
  na.omit()

#weighted mod
WMV_weighted <- readRDS("./cv_sample_A_global_vols_mods/model_objs/WMV_lambdaNULL_all_sex_mod.rds")

wp.taki(xvar=df$logAge_days, resid=WMV_weighted$residuals, n.inter=4)
summary(WMV_weighted)

WMV_unweighted <- gamlss(formula = WMV ~ pb(sexMale_x_logAge, control = pb.control(order = 3)) + pb(logAge_days, control = pb.control(order = 3)) + random(study_site), 
                         sigma.formula = ~pb(sexMale_x_logAge, control = pb.control(order = 3)) + pb(logAge_days, control = pb.control(order = 3)) + random(study_site),  
                         nu.formula = ~sexMale + fs_version_GM, family = GG, data = df, start.from = mod_list[[1]],  
                         control = gamlss.control(n.cyc = 200, nu.step = 0.25), trace = FALSE) 