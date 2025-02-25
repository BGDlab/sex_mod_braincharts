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
#penalty lambda on order 3
gamlss_3lambda <- function(pheno, lambda=NULL, fs_ver, fam="GG"){
  result <- tryCatch({
    gamlss_RSformula <-paste("gamlss(formula =", pheno, "~ pb(sexMale_x_logAge, control=pb.control(order=3, lambda=", lambda, ")) + pb(logAge_days, control=pb.control(order=3, lambda=", lambda, ")) + random(study_site) +", fs_ver, ",",
                             "sigma.formula = ~ pb(sexMale_x_logAge, control=pb.control(order=3, lambda=", lambda, ")) + pb(logAge_days, control=pb.control(order=3, lambda=", lambda, ")) + random(study_site) +", fs_ver, ",",
                             "nu.formula = ~ 1, control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste("gamlss(formula =", pheno, "~ pb(sexMale_x_logAge, control=pb.control(order=3, lambda=", lambda, ")) + pb(logAge_days, control=pb.control(order=3, lambda=", lambda, ")) + random(study_site) +", fs_ver, ",",
                               "sigma.formula = ~ pb(sexMale_x_logAge, control=pb.control(order=3, lambda=", lambda, ")) + pb(logAge_days, control=pb.control(order=3, lambda=", lambda, ")) + random(study_site) +", fs_ver, ",",
                               "nu.formula = ~ 1, method=CG(), control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
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

#gamlss_3pb
#penalized splines with order=3
gamlss_3pb <- function(pheno, fs_ver, fam="GG"){
  result <- tryCatch({
    gamlss_RSformula <-paste("gamlss(formula =", pheno, "~ pb(sexMale_x_logAge, control=pb.control(order=3)) + pb(logAge_days, control=pb.control(order=3)) + random(study_site) +", fs_ver, ",",
                             "sigma.formula = ~ pb(sexMale_x_logAge, control=pb.control(order=3)) + pb(logAge_days, control=pb.control(order=3)) + random(study_site) +", fs_ver, ",",
                             "nu.formula = ~ 1, control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste("gamlss(formula =", pheno, "~ pb(sexMale_x_logAge, control=pb.control(order=3)) + pb(logAge_days, control=pb.control(order=3)) + random(study_site) +", fs_ver, ",",
                               "sigma.formula = ~ pb(sexMale_x_logAge, control=pb.control(order=3)) + pb(logAge_days, control=pb.control(order=3)) + random(study_site) +", fs_ver, ",",
                               "nu.formula = ~ 1, method=CG(), control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
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

#gamlss_3pb_etiv
#penalized splines with order=3
gamlss_3pb_etiv <- function(pheno, fs_ver, fam="GG"){
  result <- tryCatch({
    gamlss_RSformula <-paste("gamlss(formula =", pheno, "~ pb(sexMale_x_logAge, control=pb.control(order=3)) + pb(logAge_days, control=pb.control(order=3)) + random(study_site) + eTIV +", fs_ver, ",",
                             "sigma.formula = ~ pb(sexMale_x_logAge, control=pb.control(order=3)) + pb(logAge_days, control=pb.control(order=3)) + random(study_site) + eTIV +", fs_ver, ",",
                             "nu.formula = ~ 1, control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste("gamlss(formula =", pheno, "~ pb(sexMale_x_logAge, control=pb.control(order=3)) + pb(logAge_days, control=pb.control(order=3)) + random(study_site) + eTIV +", fs_ver, ",",
                               "sigma.formula = ~ pb(sexMale_x_logAge, control=pb.control(order=3)) + pb(logAge_days, control=pb.control(order=3)) + random(study_site) + eTIV +", fs_ver, ",",
                               "nu.formula = ~ 1, method=CG(), control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
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
