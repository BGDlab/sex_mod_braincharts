#functions to fit various gamlss models

library(gamlss)
#library(gamlss2)

#gamlss_knots()
#basic_ns_knots
gamlss_knots <- function(pheno, knots=NULL, sigma_knots=NULL, fam= "GG"){
  result <- tryCatch({
    gamlss_RSformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + sexMale + fs_version + study_site,",
                             "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + sexMale + fs_version + study_site,",
                             "nu.formula = ~ 1, control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste("gamlss(formula =", pheno, "~ ns(logAge_days, knots = c(", knots, ")) + sexMale + fs_version + study_site,",
                               "sigma.formula = ~ ns(logAge_days, knots = c(", sigma_knots, ")) + sexMale + fs_version + study_site,",
                               "nu.formula = ~ 1, method=CG(), control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
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

#gamlss_3lambda
#penalty lambda on order 3 - using with default lambda=NULL -> model will select
gamlss_3lambda <- function(pheno, lambda=NULL, fs_ver, fs_moment=c("both", "mu", "none"), fam="GG"){
  
  fs_moment <- match.arg(fs_moment)
  
  #define formulas for each moment
  
  mu_base <- paste("gamlss(formula =", pheno, "~ pb(sexMale_x_logAge, lambda=", lambda, ", control=pb.control(order=3)) + pb(logAge_days, lambda=", lambda, ", control=pb.control(order=3)) + random(study_site)")
  
  if (fs_moment == "both" | fs_moment == "mu"){
    mu_form <- paste(mu_base, "+", fs_ver, ",") #add fs_version term if needed
  } else {
    mu_form <- paste(mu_base, ",") #or just comma
  }
  
  sig_base <- paste("sigma.formula = ~ pb(sexMale_x_logAge, lambda=", lambda, ", control=pb.control(order=3)) + pb(logAge_days, lambda=", lambda, ", control=pb.control(order=3)) + random(study_site)")
  
  if (fs_moment == "both") {
    sig_form <- paste(sig_base, "+", fs_ver, ",")
  } else {
    sig_form <- paste(",")
  }
  
  nu_form <- paste("nu.formula = ~ 1,")
  control <- paste("control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
  
  #try methods
  
  result <- tryCatch({
    gamlss_RSformula <-paste(mu_form, sig_form, nu_form, control)
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste(mu_form, sig_form, nu_form, "method=CG()", control)
      eval(parse(text = gamlss_CGformula))
      
      #if CG also fails, return NULL
    }, error = function(e2) {
      message(e2$message, ", returning NULL")
      return(NULL)
    })
  } , finally = {
    message("done")
  } )
  return(result)
}


