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
gamlss_3lambda <- function(pheno, lambda=NULL, 
                           fs_ver, fs_moment=c("both", "mu", "none", "all"), 
                           fam="GG",
                           nu_form="1",
                           start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  
  #define formulas for each moment
  mu_base <- paste(
    "gamlss(formula =", pheno, "~",
    make_pb("sexMale_x_logAge", lambda), "+",
    make_pb("logAge_days", lambda), "+ random(study_site)"
  )
  
  if (fs_moment == "both" | fs_moment == "mu"){
    mu_form <- paste(mu_base, "+", fs_ver) #add fs_version term if needed
  } else {
    mu_form <- paste(mu_base) #or just comma
  }
  
  sig_base <- paste(
    "sigma.formula = ~",
    make_pb("sexMale_x_logAge", lambda), "+",
    make_pb("logAge_days", lambda), "+ random(study_site)"
  )
  
  if (fs_moment == "both") {
    sig_form <- paste(sig_base, "+", fs_ver)
  } else {
    sig_form <- sig_base
  }
  
  if (fs_moment == "all") {
    nu_form <- paste("nu.formula = ~", nu_form," + ", fs_ver)
  } else {
    nu_form <- paste("nu.formula = ~", nu_form)
  }
  
  if (is.null(start.from)) {
    control <- paste("control = gamlss.control(n.cyc = 200), family =", fam, ", data= df, trace = FALSE)")
  } else if (is.gamlss(start.from)) {
    control <- paste0("start.from = ", start.from,
                      ", control = gamlss.control(n.cyc = 200), family =", fam,
                      ", data= df, trace = FALSE)")
  } else {
    stop("start.from arg must be gamlss model")
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
  return(result)
}

#gamlss_3lambda_rep
#rep model with lambdas from another model obj
gamlss_3lambda_rep <- function(og_mod, null_mod=TRUE){
  
  pheno <- paste(og_mod$mu.formula)[2]
  fam <- og_mod$family[1]

  #define formulas for each moment
  
  #MU
  mu_base <- paste(og_mod$mu.formula)[3]
  mu_lambdas <- og_mod$mu.lambda
  
  #update lambdas
  if (null_mod == FALSE){
    mu_base <- sub("sexMale_x_logAge, lambda = *,", paste("sexMale_x_logAge, lambda =", mu_lambdas[1], ","), mu_base)
  } else if (null_mod == TRUE) {
    mu_base <- sub("pb\\(sexMale_x_logAge,\\s*lambda\\s*=\\s*[^,]+,\\s*control\\s*=\\s*pb\\.control\\(order\\s*=\\s*3\\)\\)\\s*\\+\\s*", 
                   "", mu_base)
  }
  
  mu_base <- sub("logAge_days, lambda = *,", paste("logAge_days, lambda =", mu_lambdas[2], ","), mu_base)
  mu_base <- sub("random(study_site", paste("random(study_site, lambda =", mu_lambdas[3], ")"), mu_base)
  
  mu_form <- paste("gamlss(formula =", pheno, "~", mu_base, ",")
  
  #SIGMA
  sig_base <- paste(og_mod$sigma.formula)[2]
  sig_lambdas <- og_mod$sigma.lambda
  
  #update lambdas
  if (null_mod == FALSE){
    sig_base <- sub("sexMale_x_logAge, lambda = *,", paste("sexMale_x_logAge, lambda =", sig_lambdas[1], ","), sig_base)
  } else if (null_mod == TRUE) {
    sig_base <- sub("pb\\(sexMale_x_logAge,\\s*lambda\\s*=\\s*[^,]+,\\s*control\\s*=\\s*pb\\.control\\(order\\s*=\\s*3\\)\\)\\s*\\+\\s*", 
                   "", sig_base)
  }
  sig_base <- sub("logAge_days, lambda = *,", paste("logAge_days, lambda =", sig_lambdas[2], ","), sig_base)
  sig_base <- sub("random(study_site", paste("random(study_site, lambda =", sig_lambdas[3], ")"), sig_base)
  
  sig_form <- paste("sigma.formula = ~", sig_base, ",")
  
  #NU
  nu_base <- paste(og_mod$sigma.formula)[2]
  nu_form <- paste("nu.formula = ~", nu_base, ",")
  
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

#helper fun
make_pb <- function(var, lambda) {
  if (is.null(lambda)) {
    paste0("pb(", var, ", control=pb.control(order=3))")
  } else {
    paste0("pb(", var, ", lambda=", lambda, ", control=pb.control(order=3))")
  }
}
