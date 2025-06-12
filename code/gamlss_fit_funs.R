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
                           weight= FALSE,
                           nu_form="1",
                           start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  
  #define formulas for each moment
  mu_base <- paste(
    "safe_gamlss(formula =", pheno, "~",
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

#gamlss_3lambda_rep
#rep model with lambdas from another model obj
gamlss_3lambda_rep <- function(og_mod, 
                               null_mod=FALSE, 
                               start.from=NULL,
                               weight=FALSE){
  
  pheno <- paste0(og_mod$mu.formula)[[2]]
  fam <- og_mod$family[1]

  #define formulas for each moment
  
  #MU
  mu_base <- paste0(og_mod$mu.formula)[[3]]
  mu_lambdas <- og_mod$mu.lambda

  #update lambdas
  if (null_mod == FALSE){
    mu_base <- sub('sexMale_x_logAge, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                   paste0("sexMale_x_logAge, lambda =", mu_lambdas[1]), mu_base)
  } else if (null_mod == TRUE) {
    mu_base <- rm_sexage(mu_base)
  }
  
  mu_base <- sub('logAge_days, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                 paste0("logAge_days, lambda =", mu_lambdas[2]), mu_base)
  mu_base <- sub("random\\(study_site\\)", paste0("random(study_site, lambda =", mu_lambdas[3], ")"), mu_base)
  
  mu_form <- paste0("safe_gamlss(formula =", pheno, "~", mu_base)
  
  #SIGMA
  sig_base <- paste0(og_mod$sigma.formula)[[2]]
  sig_lambdas <- og_mod$sigma.lambda
  
  #update lambdas
  if (null_mod == FALSE){
    sig_base <- sub('sexMale_x_logAge, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                    paste0("sexMale_x_logAge, lambda =", sig_lambdas[1]), sig_base)
  } else if (null_mod == TRUE) {
    sig_base <- rm_sexage(sig_base)
  }

  sig_base <- sub('logAge_days, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                  paste0("logAge_days, lambda =", sig_lambdas[2]), sig_base)
  sig_base <- sub("random\\(study_site\\)", paste0("random(study_site, lambda =", sig_lambdas[3], ")"), sig_base)
  
  sig_form <- paste0("sigma.formula = ~", sig_base)
  
  #NU
  nu_base <- paste0(og_mod$nu.formula)[[2]]
  if (!is.null(og_mod$nu.lambda)){
    nu_lambdas <- og_mod$nu.lambda
    nu_base <- sub('logAge_days, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                    paste0("logAge_days, lambda =", nu_lambdas[1]), nu_base)
  }
  nu_form <- paste0("nu.formula = ~", nu_base)
  
  control <- paste("control = gamlss.control(n.cyc = 200, nu.step=0.25), family =", og_mod$family[[1]], ", data= df, trace = FALSE)")

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

#helper fun
make_pb <- function(var, lambda) {
  if (is.null(lambda)) {
    paste0("pb(", var, ", method='GAIC', k=log(nrow(df)), control=pb.control(order=3))")
  } else {
    paste0("pb(", var, ", lambda=", lambda, ", method='GAIC', k=log(nrow(df)), control=pb.control(order=3))")
  }
}

#gamlss_3lambda_etiv
#penalty lambda on order 3 - using with default lambda=NULL -> model will select
gamlss_3lambda_etiv <- function(pheno, lambda=NULL,
                                total_var, total_moment=c("both", "mu", "none", "all"),
                                fs_ver, fs_moment=c("both", "mu", "none", "all"), 
                                fam="GG", 
                                weight= FALSE,
                                nu_form="1",
                                start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  total_moment <- match.arg(total_moment)
  
  #define formulas for each moment
  #MU BASE
  mu_base <- paste(
    "safe_gamlss(formula =", pheno, "~",
    make_pb("sexMale_x_logAge", lambda), "+",
    make_pb("logAge_days", lambda), "+ sexMale + random(study_site)"
  )
  
  #add fs_version
  if (fs_moment != "none"){
    mu_base <- paste(mu_base, "+", fs_ver)
  }
  
  #add etiv/total pheno
  if (total_moment != "none"){
    mu_base <- paste(mu_base, "+", make_pb(total_var, lambda))
  }
  
  #SIGMA BASE
  sig_base <- paste(
    "sigma.formula = ~",
    make_pb("sexMale_x_logAge", lambda), "+",
    make_pb("logAge_days", lambda), "+ sexMale + random(study_site)"
  )
  
  if (fs_moment == "both" | fs_moment == "all") {
    sig_base <- paste(sig_base, "+", fs_ver)
  }
  
  if (total_moment == "both" | total_moment == "all") {
    sig_base <- paste(sig_base, "+", make_pb(total_var, lambda))
  }
  
  #NU BASE
  nu_base <- paste("nu.formula = ~", nu_form)
  
  if (fs_moment == "all") {
    nu_base <- paste(nu_base," + ", fs_ver)
  } 
  
  if (total_moment == "all") {
    nu_base <- paste(nu_base," + ", total_var)
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
    gamlss_RSformula <-paste(mu_base, sig_base, nu_base, control, sep=", ")
    print(gamlss_RSformula)
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste(mu_base, sig_base, nu_base, "method=CG()", control, sep=", ")
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
    control <- sub("nu\\.step\\s*=\\s*[-+]?[0-9]*\\.?[0-9]+", "nu.step = 0.00000000001", control)
    
    result <- tryCatch({
      gamlss_RSformula <-paste(mu_base, sig_base, nu_base, control, sep=", ")
      print(gamlss_RSformula)
      
      eval(parse(text = gamlss_RSformula))
      
    } , warning = function(w) {
      message("warning")
      eval(parse(text = gamlss_RSformula))
      
    } , error = function(e) {
      message(e$message, ", trying method=CG()")
      tryCatch({
        gamlss_CGformula <-paste(mu_base, sig_base, nu_base, "method=CG()", control, sep=", ")
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

#gamlss_4param
gamlss_4param <- function(pheno, 
                          mu_base,
                          sig_base,
                          nu_base,
                          tau_base,
                          fam="BCT",
                          start.from=NULL){
  

  #define formulas for each moment
  mu_form <- paste("gamlss(formula =", pheno, "~", mu_base)
  sig_form <- paste("sigma.formula = ~", sig_base)
  nu_form <- paste("nu.formula = ~", nu_base)
  tau_form <- paste("ta.formula = ~", tau_base)
  
  
  control <- paste("control = gamlss.control(n.cyc = 200, nu.step=0.25), family =", fam, ", data= df, trace = FALSE)")
  
  if (!is.null(start.from)) {
    control <- paste0("start.from = ", start.from,", ", control)
  }
  
  #try methods
  
  result <- tryCatch({
    gamlss_RSformula <-paste(mu_form, sig_form, nu_form, tau_form, control, sep=", ")
    print(gamlss_RSformula)
    
    eval(parse(text = gamlss_RSformula))
    
  } , warning = function(w) {
    message("warning")
    eval(parse(text = gamlss_RSformula))
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste(mu_form, sig_form, nu_form, tau_form, "method=CG()", control, sep=", ")
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
    control <- sub("nu\\.step\\s*=\\s*[-+]?[0-9]*\\.?[0-9]+", "nu.step = 0.00000000001", control)
    control <- sub("tau\\.step\\s*=\\s*[-+]?[0-9]*\\.?[0-9]+", "tau.step = 0.00000000001", control)
    
    result <- tryCatch({
      gamlss_RSformula <-paste(mu_form, sig_form, nu_form, tau_form, control, sep=", ")
      print(gamlss_RSformula)
      
      eval(parse(text = gamlss_RSformula))
      
    } , warning = function(w) {
      message("warning")
      eval(parse(text = gamlss_RSformula))
      
    } , error = function(e) {
      message(e$message, ", trying method=CG()")
      tryCatch({
        gamlss_CGformula <-paste(mu_form, sig_form, nu_form, tau_form, "method=CG()", control, sep=", ")
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


#written with help from GPT, trying to throw error instead of silently returning NULL model
safe_gamlss <- function(...) {
  mod <- gamlss(...)
  
  # Check for NULL coefficients
  null_mu <- is.null(coef(mod, what = "mu"))
  null_sigma <- is.null(coef(mod, what = "sigma"))
  null_nu <- is.null(coef(mod, what = "nu"))
  
  if (null_mu && null_sigma && null_nu) {
    stop("Model fit failed: coefficients are NULL")
  }
  
  return(mod)
}

#try not log-scaling age
gamlss_age <- function(pheno, lambda=NULL, 
                           fs_ver, fs_moment=c("both", "mu", "none", "all"), 
                           fam="GG",
                           weight= FALSE,
                           nu_form="1",
                           start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  
  #define formulas for each moment
  mu_base <- paste(
    "safe_gamlss(formula =", pheno, "~",
    make_pb("sexMale_x_age", lambda), "+",
    make_pb("age_days", lambda), "+ sexMale + random(study_site)"
  )
  
  if (fs_moment != "none"){
    mu_form <- paste(mu_base, "+", fs_ver) #add fs_version term if needed
  } else {
    mu_form <- paste(mu_base) #or just comma
  }
  
  sig_base <- paste(
    "sigma.formula = ~",
    make_pb("sexMale_x_age", lambda), "+",
    make_pb("age_days", lambda), "+ sexMale + random(study_site)"
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

#cs() to smooth instead of pb()
gamlss_cs <- function(pheno,
                      fs_ver, 
                      fs_moment=c("both", "mu", "none", "all"), 
                      fam="GG",
                      weight= FALSE,
                      nu_form="1",
                      start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  
  #define formulas for each moment
  mu_base <- paste(
    "safe_gamlss(formula =", pheno, "~ cs(sexMale_x_logAge) + cs(logAge_days)+ sexMale + random(study_site)"
  )
  
  if (fs_moment != "none"){
    mu_form <- paste(mu_base, "+", fs_ver) #add fs_version term if needed
  } else {
    mu_form <- paste(mu_base) #or just comma
  }
  
  sig_base <- paste(
    "sigma.formula = ~ cs(sexMale_x_logAge) + cs(logAge_days) + sexMale + random(study_site)"
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

#cs() to smooth with unscaled age
gamlss_csage <- function(pheno,
                      fs_ver, 
                      fs_moment=c("both", "mu", "none", "all"), 
                      fam="GG",
                      weight= FALSE,
                      nu_form="1",
                      start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  
  #define formulas for each moment
  mu_base <- paste(
    "safe_gamlss(formula =", pheno, "~ cs(sexMale_x_age) + cs(age_days)+ sexMale + random(study_site)"
  )
  
  if (fs_moment != "none"){
    mu_form <- paste(mu_base, "+", fs_ver) #add fs_version term if needed
  } else {
    mu_form <- paste(mu_base) #or just comma
  }
  
  sig_base <- paste(
    "sigma.formula = ~ cs(sexMale_x_age) + cs(age_days) + sexMale + random(study_site)"
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


#rm
rm_sexage <- function(formula_string) {
  sub(
    pattern = 'pb\\(sexMale_x_logAge, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\), control = pb\\.control\\(order = 3\\)\\) \\+\\s*',
    replacement = '',
    x = formula_string
  )
}

#scaling/unscaling phenotype functions
log_scale <- function(x){log(x + 5, base=10)}
unscale <- function(x){10^x - 5}
