#quick test for adding 1000 to phenos

set.seed(9999)

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

source("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/gamlss_fit_funs.R")

df_full <- fread("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/data/cv_sample_A.csv", stringsAsFactors=TRUE)
pheno_list <- readRDS("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/pheno_lists/global_vols.rds")

gamlss_3lambda1k <- function(pheno, 
                             lambda=NULL, 
                             fs_ver, fs_moment=c("both", "mu", "none", "all"), 
                             fam="BCCG",
                             weight= FALSE,
                             nu_form="1",
                             start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  
  #define formulas for each moment
  mu_base <- paste(
    "safe_gamlss(formula =", pheno, "~ 1000 +",
    make_pb("sexMale_x_logAge", lambda), "+",
    make_pb("logAge_days", lambda), "+ sexMale + random(study_site)"
  )
  
  if (fs_moment != "none"){
    mu_form <- paste(mu_base, "+", fs_ver) #add fs_version term if needed
  } else {
    mu_form <- paste(mu_base) #or just comma
  }
  
  sig_base <- paste(
    "sigma.formula = ~",
    make_pb("sexMale_x_logAge", lambda), "+",
    make_pb("logAge_days", lambda), "+ sexMale + random(study_site)"
  )
  
  if (fs_moment == "both" | fs_moment == "all") {
    sig_form <- paste(sig_base, "+", fs_ver)
  } else {
    sig_form <- sig_base
  }
  
  if (fs_moment == "all") {
    nu_form <- paste("nu.formula = ~", nu_form," + ", fs_ver)
  } else {
    nu_form <- paste("nu.formula = ~", nu_form)
  }
  
  control <- paste("control = gamlss.control(n.cyc = 200, nu.step=0.25), family =", fam, ", data= df, trace = FALSE)")
  
  if (!is.null(start.from)) {
    control <- paste0("start.from = ", start.from,", ", control)
  }
  
  if (weight==TRUE) {
    control <- paste0("weights = weight, ", control)
  }
  
  #try methods
  
  result <- tryCatch({
    gamlss_RSformula <-paste(mu_form, sig_form, nu_form, control, sep=", ")
    print(gamlss_RSformula)
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste(mu_form, sig_form, nu_form, "method=CG()", control, sep=", ")
      eval(parse(text = gamlss_CGformula))
      
      #if CG also fails, return NULL
    }, error = function(e2) {
      message(e2$message, ", returning NULL")
      return(NULL)
    })
  } , finally = {
    message("done")
  } )
  
  #if needed, try again with tiny nu.step
  if(is.null(result)){
    control <- sub("mu\\.step\\s*=\\s*[-+]?[0-9]*\\.?[0-9]+", "mu.step = 0.00000000001", control)
    control <- sub("sigma\\.step\\s*=\\s*[-+]?[0-9]*\\.?[0-9]+", "sigma.step = 0.00000000001", control)
    control <- sub("nu\\.step\\s*=\\s*[-+]?[0-9]*\\.?[0-9]+", "nu.step = 0.00000000001", control)
    
    result <- tryCatch({
      gamlss_RSformula <-paste(mu_form, sig_form, nu_form, control, sep=", ")
      print(gamlss_RSformula)
      
      eval(parse(text = gamlss_RSformula))
      
    } , warning = function(w) {
      message("warning")
      eval(parse(text = gamlss_RSformula))
      
    } , error = function(e) {
      message(e$message, ", trying method=CG()")
      tryCatch({
        gamlss_CGformula <-paste(mu_form, sig_form, nu_form, "method=CG()", control, sep=", ")
        eval(parse(text = gamlss_CGformula))
        
        #if CG also fails, return NULL
      }, error = function(e2) {
        message(e2$message, ", returning NULL")
        return(NULL)
      })
    } , finally = {
      message("done")
    } )
    
  }
  
  return(result)
}

for (pheno in pheno_list){
  pheno <- as.character(pheno)
  print(pheno)
  
  pheno_sym <- sym(pheno)
  
  df <- df_full %>%
    dplyr::select(all_of(c(pheno, "fs_version_GM", "logAge_days", "sexMale_x_logAge", "sexMale", "study_site"))) %>%
    mutate(!!pheno_sym := (!!pheno_sym + 1000)) %>% #transform    
    na.omit() %>%
    trunc_coverage("logAge_days")
    
  test_mod <- gamlss_3lambda1k(pheno,
                               lambda=NULL, 
                               fs_ver = "fs_version_GM",
                               fs_moment="both",
                               fam="BCCG", 
                               weight= FALSE,
                               nu_form="logAge_days",
                               start.from=NULL)
  
  saveRDS(test_mod, paste0("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/", pheno, "_test_mod.rds"))
  
  fan_plot <- make_centile_fan(test_mod, df, x_var="logAge_days", color_var="sexMale")
  ggsave(file=paste0("/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/", pheno, "_test_centiles.png"), fan_plot)
  
}