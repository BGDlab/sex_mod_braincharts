#functions to fit various gamlss models

library(gamlss)
#library(gamlss2)

#gamlss_lambda
gamlss_lambda <- function(pheno, lambda=NULL, 
                           fs_ver, fs_moment=c("both", "mu", "none", "all"), 
                           fam="GG",
                           weight=NULL,
                           nu_form="1",
                           start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  
  #define formulas for each moment
  mu_base <- paste(
    "safe_gamlss_old(formula =", pheno, "~",
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
  
  control <- paste("control = gamlss.control(n.cyc=400, nu.step=0.25), family =", fam, ", data= df, trace = FALSE)")
  
  if (!is.null(start.from)) {
    control <- paste0("start.from = ", start.from,", ", control)
  }
  
  if (!is.null(weight)) {
    control <- paste0("weights = df$", weight, ",", control)
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

#gamlss_lambda_rep
#rep model with lambdas from another model obj
gamlss_lambda_rep <- function(og_mod, 
                               null_mod=c("false", "true", "allSex"),
                               keep_lambdas=TRUE,
                               start.from=NULL,
                               weight=NULL,
                               n.cyc=400){
  
  pheno <- paste0(og_mod$mu.formula)[[2]]
  fam <- og_mod$family[1]
  null_mod <- match.arg(null_mod)

  #define formulas for each moment
  
  #MU
  mu_base <- paste0(og_mod$mu.formula)[[3]]
  mu_lambdas <- og_mod$mu.lambda

  #update lambdas
  if (keep_lambdas == TRUE & null_mod == "false"){
    mu_base <- sub('sexMale_x_logAge, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                   paste0("sexMale_x_logAge, lambda =", mu_lambdas[1]), mu_base)
  } else if (null_mod != "false") {
    mu_base <- rm_sexage(mu_base)
    if (null_mod == "allSex"){
      #remove sex intercept
      mu_base <- sub('+ sexMale', '', mu_base)
    }
  } 
  
  if (keep_lambdas == TRUE){
    mu_base <- sub('logAge_days, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                   paste0("logAge_days, lambda =", mu_lambdas[2]), mu_base)
    mu_base <- sub("random\\(study_site\\)", paste0("random(study_site, lambda =", mu_lambdas[3], ")"), mu_base)
    
    #replace TBV or other pbs if needed
    if (length(mu_lambdas)==4){
      mu_base <- sub('method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                     paste0("lambda =", mu_lambdas[4]), mu_base)
    }
  }
  
  mu_form <- paste0("safe_gamlss_old(formula =", pheno, "~", mu_base)
  
  #SIGMA
  sig_base <- paste0(og_mod$sigma.formula)[[2]]
  sig_lambdas <- og_mod$sigma.lambda
  
  #update lambdas
  if (keep_lambdas == TRUE & null_mod == "false"){
    sig_base <- sub('sexMale_x_logAge, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                    paste0("sexMale_x_logAge, lambda =", sig_lambdas[1]), sig_base)
  } else if (null_mod != "false") {
    sig_base <- rm_sexage(sig_base)
    if (null_mod == "allSex"){
      #remove sex intercept
      sig_base <- sub('+ sexMale', '', sig_base)
    }
  }
  
  if (keep_lambdas==TRUE){
    sig_base <- sub('logAge_days, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                    paste0("logAge_days, lambda =", sig_lambdas[2]), sig_base)
    sig_base <- sub("random\\(study_site\\)", paste0("random(study_site, lambda =", sig_lambdas[3], ")"), sig_base)
    
    #replace TBV or other pbs if needed
    if (length(sig_lambdas)==4){
      sig_base <- sub('method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                     paste0("lambda =", sig_lambdas[4]), sig_base)
    }
  }
  
  sig_form <- paste0("sigma.formula = ~", sig_base)
  
  #NU
  nu_base <- paste0(og_mod$nu.formula)[[2]]
  if (!is.null(og_mod$nu.lambda) & keep_lambdas==TRUE){
    nu_lambdas <- og_mod$nu.lambda
    nu_base <- sub('logAge_days, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)', 
                    paste0("logAge_days, lambda =", nu_lambdas[1]), nu_base)
  }
  
  if (null_mod == "allSex"){
    #remove sex intercept
    nu_base <- sub('+ sexMale', '', nu_base)
  }
  
  nu_form <- paste0("nu.formula = ~", nu_base)
  
  control <- paste("control = gamlss.control(n.cyc=", n.cyc,", nu.step=0.25), family =", og_mod$family[[1]], ", data= df, trace = FALSE)")

  if (!is.null(start.from)) {
    control <- paste0("start.from = ", start.from,", ", control)
  }
  
  if (!is.null(weight)) {
    control <- paste0("weights = df$", weight, ",", control)
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
  
  #if needed, try one last time with mixed method
  if(is.null(result)){
    result <- tryCatch({
      gamlss_Mformula <-paste(mu_form, sig_form, nu_form, "method=mixed(10,500)", control, sep=", ")
      print(gamlss_Mformula)
      
      eval(parse(text = gamlss_Mformula))
      
    } , warning = function(w) {
      message("warning")
      eval(parse(text = gamlss_Mformula))
      
    } , error = function(e) {
      message(e$message, ", trying method=CG()")
      tryCatch({
        gamlss_Mformula <-paste(mu_form, sig_form, nu_form, "method=mixed(10,800)", control, sep=", ")
        eval(parse(text = gamlss_Mformula))
        
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
    paste0("pb(", var, ", method='GAIC', k=log(nrow(df)))")
  } else {
    paste0("pb(", var, ", lambda=", lambda, ", method='GAIC', k=log(nrow(df)))")
  }
}

#gamlss_lambda_etiv
#using with default lambda=NULL -> model will select
gamlss_lambda_etiv <- function(pheno, lambda=NULL,
                                total_var, total_moment=c("both", "mu", "none", "all"),
                                fs_ver, fs_moment=c("both", "mu", "none", "all"), 
                                fam="GG", 
                                weight= NULL,
                                nu_form="1",
                                start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  total_moment <- match.arg(total_moment)
  
  #define formulas for each moment
  #MU BASE
  mu_base <- paste(
    "safe_gamlss_old(formula =", pheno, "~",
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
  
  control <- paste("control = gamlss.control(n.cyc=400, nu.step=0.25), family =", fam, ", data= df, trace = FALSE)")
  
  if (!is.null(start.from)) {
    control <- paste0("start.from = ", start.from,", ", control)
  }
  
  if (!is.null(weight)) {
    control <- paste0("weights = df$", weight, ",", control)
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


#try not log-scaling age
gamlss_age <- function(pheno, lambda=NULL, 
                           fs_ver, fs_moment=c("both", "mu", "none", "all"), 
                           fam="GG",
                           weight= NULL,
                           nu_form="1",
                           start.from=NULL){
  
  fs_moment <- match.arg(fs_moment)
  
  #define formulas for each moment
  mu_base <- paste(
    "safe_gamlss_old(formula =", pheno, "~",
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
  
  control <- paste("control = gamlss.control(n.cyc=400, nu.step=0.25), family =", fam, ", data= df, trace = FALSE)")
  
  if (!is.null(start.from)) {
    control <- paste0("start.from = ", start.from,", ", control)
  }
  
  if (!is.null(weight)) {
    control <- paste0("weights = df$", weight, ",", control)
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
    pattern = 'pb\\(sexMale_x_logAge, method = \\"GAIC\\", k = log\\(nrow\\(df\\)\\)\\) \\+\\s*',
    replacement = '',
    x = formula_string
  )
}

#scaling/unscaling phenotype functions
log_scale <- function(x){log(x + 5, base=10)}
unscale <- function(x){10^x - 5}

rm_lambdas <- function(formula_string){
  text_clean <- gsub(
    pattern = "lambda\\s*=\\s*[0-9.]+\\s*(?=[),])",
    replacement = "lambda = NULL",
    x = formula_string,
    perl = TRUE
  )
  return(text_clean)
}

#tmp rollback to fix CG() and mixed() methods
safe_gamlss_old <- function(...) {
  warn_msg <- NULL
  
  
  mod <- withCallingHandlers({
    gamlss(...)
  }, warning = function(w) {
    # Capture the warning message
    warn_msg <<- w$message
    
    # Example condition: promote warnings containing "Error" or convergence issues
    if (grepl("Error", w$message, ignore.case = TRUE) ||
        grepl("converge", w$message, ignore.case = TRUE)) {
      # Turn this warning into an error
      stop(simpleError(w$message))
    }
  },
  error = function(e) {
    stop(e)  # propagate any real errors
  }
  )
  
  # Check for NULL coefficients
  null_mu <- is.null(coef(mod, what = "mu"))
  null_sigma <- is.null(coef(mod, what = "sigma"))
  
  if (null_mu && null_sigma) {
    stop("Model fit failed: coefficients are NULL")
  }
  
  #backup check
  if (mod$converged==FALSE) {
    stop("Model did not converge:", warn_msg)
  }
  
  return(mod)
}