#functions to fit various gamlss models

library(gamlss)
library(gamlss2)

#gamlss_knots()
#basic_ns_knots
gamlss_knots <- function(pheno, knots=NULL, sigma_knots=NULL){
  result <- tryCatch({
    gamlss_RSformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + sexMale + fs_version + study_site,",
                             "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + sexMale + fs_version + study_site,",
                             "nu.formula = ~ 1, control = gamlss.control(n.cyc = 200), family = GG, data= df, trace = FALSE)")
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + sexMale + fs_version + study_site,",
                               "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + sexMale + fs_version + study_site,",
                               "nu.formula = ~ 1, method=CG(), control = gamlss.control(n.cyc = 200), family = GG, data= df, trace = FALSE)")
      eval(parse(text = gamlss_CGformula))
      
      #if CG also fails, return NULL
    }, error = function(e2) {
      message("second error, returning NULL")
      return(NULL)
    })
  } , finally = {
    message("done")
  } )
  return(result)
}

#gamlss_mod_knots()
#same as gamlss_knots() but with age*sex effects
gamlss_mod_knots <- function(pheno, knots=NULL, sigma_knots=NULL){
  result <- tryCatch({
    gamlss_RSformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + ns(sexMale_x_logAge, knots = c(", knots, ")) + sexMale + fs_version + study_site,",
                             "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + ns(sexMale_x_logAge, knots = c(", sigma_knots, ")) + sexMale + fs_version + study_site,",
                             "nu.formula = ~ 1, control = gamlss.control(n.cyc = 200), family = GG, data= df, trace = FALSE)")
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + ns(sexMale_x_logAge, knots = c(", knots, ")) + sexMale + fs_version + study_site,",
                               "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + ns(sexMale_x_logAge, knots = c(", sigma_knots, ")) + sexMale + fs_version + study_site,",
                               "nu.formula = ~ 1, method=CG(), control = gamlss.control(n.cyc = 200), family = GG, data= df, trace = FALSE)")
      eval(parse(text = gamlss_CGformula))
      
      #if CG alos fails, return NULL
    }, error = function(e2) {
      message(e2$message, ", returning NULL")
      return(NULL)
    })
  } , finally = {
    message("done")
  } )
  return(result)
}

#gamlss_mod_nofs()
#same as gamlss_mod_knots() but without fs_version term
gamlss_mod_nofs <- function(pheno, knots=NULL, sigma_knots=NULL){
  result <- tryCatch({
    gamlss_RSformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + ns(sexMale_x_logAge, knots = c(", knots, ")) + sexMale + study_site,",
                             "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + ns(sexMale_x_logAge, knots = c(", sigma_knots, ")) + sexMale + study_site,",
                             "nu.formula = ~ 1, control = gamlss.control(n.cyc = 200), family = GG, data= df, trace = FALSE)")
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + ns(sexMale_x_logAge, knots = c(", knots, ")) + sexMale + study_site,",
                               "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + ns(sexMale_x_logAge, knots = c(", sigma_knots, ")) + sexMale + study_site,",
                               "nu.formula = ~ 1, method=CG(), control = gamlss.control(n.cyc = 200), family = GG, data= df, trace = FALSE)")
      eval(parse(text = gamlss_CGformula))
      
      #if CG alos fails, return NULL
    }, error = function(e2) {
      message(e2$message, ", returning NULL")
      return(NULL)
    })
  } , finally = {
    message("done")
  } )
  return(result)
}

#gamlss2_mod_knots()
#same as gamlss_mod_knots() but using gamlss2 - NOT READY YET
# gamlss2_mod_knots <- function(pheno, knots=NULL, sigma_knots=NULL){
#   result <- tryCatch({
#     gamlss_RSformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + ns(sexMale_x_logAge, knots = c(", knots, ")) + sexMale + fs_version + study_site,",
#                              "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + ns(sexMale_x_logAge, knots = c(", sigma_knots, ")) + sexMale + fs_version + study_site,",
#                              "nu.formula = ~ 1, control = gamlss.control(n.cyc = 200), family = GG, data= df, trace = FALSE)")
#     
#     eval(parse(text = gamlss_RSformula))
#     
#   } , warning = function(w) {
#     message("warning")
#     eval(parse(text = gamlss_RSformula))
#     
#   } , error = function(e) {
#     message(e$message, ", trying method=CG()")
#     tryCatch({
#       gamlss_CGformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + ns(sexMale_x_logAge, knots = c(", knots, ")) + sexMale + fs_version + study_site,",
#                                "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + ns(sexMale_x_logAge, knots = c(", sigma_knots, ")) + sexMale + fs_version + study_site,",
#                                "nu.formula = ~ 1, method=CG(), control = gamlss.control(n.cyc = 200), family = GG, data= df, trace = FALSE)")
#       eval(parse(text = gamlss_CGformula))
#       
#       #if CG alos fails, return NULL
#     }, error = function(e2) {
#       message(e2$message, ", returning NULL")
#       return(NULL)
#     })
#   } , finally = {
#     message("done")
#   } )
#   return(result)
# }