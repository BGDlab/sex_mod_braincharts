#functions to fit various gamlss models

library(gamlss)
#library(gamlss2)

### SEX CODING (sex_pvc branch) ###
#sex x age term in mu & sigma: pb(logAge_days) + pvc(logAge_days, by = sexEffect), sexEffect = -1 female / +1 male
#(in place of pb(sexMale_x_logAge) + pb(logAge_days)).
#- pvc() centres a numeric `by` on its mean, so the term is beta(age) * (sexEffect - mean(sexEffect)): pb(logAge_days)
#  is the curve at the sample-weighted average of the sexes and pvc the penalized sex difference (M - F = 2 * beta(age))
#- mu & sigma keep a sex main effect, as in the published models, so the null (pvc dropped) stays nested. It is
#  aliased with pvc's constant, but gamlss counts pvc as edf - 1 (its placeholder column is all zero), so the main
#  effect supplies that 1 df
#- gamlss 5.5-0 needs two fixes for a numeric `by`; see patch_pvc_numeric_by() (applied when this file is sourced)
SEX_VAR <- "sexEffect"
SEX_SPECIAL_TERM <- "sexEffect = 2 * sexMale - 1" #for sim_data(special_term =)

#add sexEffect (-1 F, +1 M) from sexMale (0/1, numeric or factor)
add_sex_effect <- function(df) {
  sex_male <- as.numeric(as.character(df$sexMale))
  stopifnot("sexMale must be 0/1" = all(sex_male %in% c(0, 1)))
  df[[SEX_VAR]] <- 2 * sex_male - 1
  df
}

#sex x age smooth terms for mu & sigma (pb before pvc: gamlss backfits smoothers in formula order)
make_sexage <- function(lambda = NULL) {
  paste(make_pb("logAge_days", lambda), "+", make_pvc(lambda), "+", SEX_VAR)
}

### MAIN FIT FUNS ###

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
    make_sexage(lambda), "+ random(study_site)"
  )

  if (fs_moment != "none"){
    mu_form <- paste(mu_base, "+", fs_ver) #add fs_version term if needed
  } else {
    mu_form <- paste(mu_base) #or just comma
  }

  sig_base <- paste(
    "sigma.formula = ~",
    make_sexage(lambda), "+ random(study_site)"
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
#each moment's formula is rebuilt from og_mod's terms: with keep_lambdas, every smoother (pb/pvc/random) gets its
#fitted lambda (gamlss stores <moment>.lambda in formula order); null_mod "true" drops the pvc sex x age term,
#"allSex" also drops the sex main effect
gamlss_lambda_rep <- function(og_mod,
                               null_mod=c("false", "true", "allSex"),
                               keep_lambdas=TRUE,
                               start.from=NULL,
                               weight=NULL,
                               n.cyc=800){

  pheno <- deparse1(og_mod$mu.formula[[2]])
  null_mod <- match.arg(null_mod)

  mu_rhs <- rep_rhs(og_mod$mu.formula, og_mod$mu.lambda, keep_lambdas, null_mod)
  sig_rhs <- rep_rhs(og_mod$sigma.formula, og_mod$sigma.lambda, keep_lambdas, null_mod)
  nu_rhs <- rep_rhs(og_mod$nu.formula, og_mod$nu.lambda, keep_lambdas, null_mod)

  if (is.null(weight)){
    mu_form <- paste0("safe_gamlss_old(formula = ", pheno, " ~ ", mu_rhs)
  } else {
    #bypass safe_gamlss when fitting with weights
    mu_form <- paste0("gamlss(formula = ", pheno, " ~ ", mu_rhs)
  }
  sig_form <- paste0("sigma.formula = ~ ", sig_rhs)
  nu_form <- paste0("nu.formula = ~ ", nu_rhs)

  control <- paste("control = gamlss.control(n.cyc=", n.cyc,", nu.step=0.25), family =", og_mod$family[[1]], ", data= df, trace = FALSE)")

  if (!is.null(start.from)) {
    control <- paste0("start.from = ", start.from,", ", control)
  }

  if (!is.null(weight)) {
    control <- paste0("weights = df$", weight, ", ", control)
  }

  #try methods

  result <- try_gamlss_methods(mu_form, sig_form, nu_form, control)

  return(result)
}

#rep_rhs: right-hand side of one moment's formula for gamlss_lambda_rep()
rep_rhs <- function(f, lambdas, keep_lambdas, null_mod) {
  tt <- terms(f)
  stopifnot(attr(tt, "intercept") == 1, is.null(attr(tt, "offset")))
  labs <- attr(tt, "term.labels")
  is_smo <- grepl("^(pb|pvc|random)\\(", labs)
  if (keep_lambdas && any(is_smo)) {
    stopifnot("one lambda per smoother" = sum(is_smo) == length(lambdas))
    labs[is_smo] <- mapply(set_lambda, labs[is_smo], lambdas, USE.NAMES = FALSE)
  }
  if (null_mod != "false") labs <- labs[!grepl("^pvc\\(", labs)]
  if (null_mod == "allSex") labs <- labs[labs != SEX_VAR]
  if (length(labs) == 0) "1" else paste(labs, collapse = " + ")
}

#set_lambda: rewrite one smoother term with a fixed lambda (dropping its lambda-selection args)
set_lambda <- function(term, lambda) {
  cl <- str2lang(term)
  cl$method <- NULL
  cl$k <- NULL
  cl$df <- NULL
  cl$lambda <- unname(lambda)
  deparse1(cl)
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
    make_sexage(lambda), "+ random(study_site)"
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
    make_sexage(lambda), "+ random(study_site)"
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

### HELPER FUNS ###

try_gamlss_methods <- function(mu_form, sig_form, nu_form, control){
  
  #first loop - RS and CG
  result <- tryCatch({
    gamlss_RSformula <-paste(mu_form, sig_form, nu_form, control, sep=", ")
    print(gamlss_RSformula)
    
    res <- eval(parse(text = gamlss_RSformula))
    
    #additional checks for when safe_gamlss() can't be used (ie fitting with weights)
    if (!isTRUE(res$converged)) {
      stop("Model did not converge")
    }
    res
    
  } , error = function(e) {
    message(e$message, ", trying method=CG()")
    tryCatch({
      gamlss_CGformula <-paste(mu_form, sig_form, nu_form, "method=CG()", control, sep=", ")
      res <- eval(parse(text = gamlss_CGformula))
      
      #additional checks for when safe_gamlss() can't be used (ie fitting with weights)
      if (!isTRUE(res$converged)) {
        stop("Model did not converge")
      }
      res
      
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
      
      res <- eval(parse(text = gamlss_RSformula))
      
      #additional checks for when safe_gamlss() can't be used (ie fitting with weights)
      if (!isTRUE(res$converged)) {
        stop("Model did not converge")
      }
      res
      
    } , error = function(e) {
      message(e$message, ", trying method=CG()")
      tryCatch({
        gamlss_CGformula <-paste(mu_form, sig_form, nu_form, "method=CG()", control, sep=", ")
        res <- eval(parse(text = gamlss_CGformula))
        
        #additional checks for when safe_gamlss() can't be used (ie fitting with weights)
        if (!isTRUE(res$converged)) {
          stop("Model did not converge")
        }
        res
        
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
        
        res <- eval(parse(text = gamlss_Mformula))
        #additional checks for when safe_gamlss() can't be used (ie fitting with weights)
        if (!isTRUE(res$converged)) {
          stop("Model did not converge")
        }
        res
        
      } , error = function(e) {
        message(e$message, ", trying more loops")
        tryCatch({
          gamlss_Mformula <-paste(mu_form, sig_form, nu_form, "method=mixed(10,800)", control, sep=", ")
          res <- eval(parse(text = gamlss_Mformula))
          
          #additional checks for when safe_gamlss() can't be used (ie fitting with weights)
          if (!isTRUE(res$converged)) {
            stop("Model did not converge")
          }
          res
          
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

#(rm_sexage() removed: gamlss_lambda_rep() now drops the pvc term via rep_rhs())

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

#helper fun
make_pb <- function(var, lambda) {
  if (is.null(lambda)) {
    paste0("pb(", var, ", method='GAIC', k=log(nrow(df)))")
  } else {
    paste0("pb(", var, ", lambda=", lambda, ", method='GAIC', k=log(nrow(df)))")
  }
}

#pvc sex x age deviation term, lambda selected by the same GAIC (k = log(n)) as the pb terms
make_pvc <- function(lambda) {
  if (is.null(lambda)) {
    paste0("pvc(logAge_days, by=", SEX_VAR, ", method='GAIC', k=log(nrow(df)))")
  } else {
    paste0("pvc(logAge_days, by=", SEX_VAR, ", lambda=", lambda, ")")
  }
}

#lambda of a model's pvc term in one moment (NA if none); <moment>.lambda is in formula order of the smoothers
pvc_lambda <- function(model, what = c("mu", "sigma")) {
  what <- match.arg(what)
  labs <- attr(terms(model[[paste0(what, ".formula")]]), "term.labels")
  smo <- labs[grepl("^(pb|pvc|random)\\(", labs)]
  i <- grep("^pvc\\(", smo)
  if (length(i) != 1) return(NA_real_)
  unname(model[[paste0(what, ".lambda")]][i])
}

#patch_pvc_numeric_by
#two gamlss 5.5-0 bugs for pvc(x, by = <numeric>), patched in this session's gamlss namespace:
#1. lambda selection: with method "GAIC", fnGAIC scores sum(w * (y - fv)^2) where fv is beta(x), fit to
#   yvar = y / by with weights w * by^2. It should score yvar, otherwise GAIC always picks the 1e7 cap (linear).
#2. predict(newdata =): pvc() centres `by` on the mean of old + new data together, so the term is refit and
#   predicted with a different centring than the fit. Re-centring on the fitting rows (the first length(y)) restores
#   the fit's centring; during fitting `by` already has mean 0, so it changes nothing there.
#stops if the source doesn't match (other gamlss version) instead of fitting with an unpatched pvc
patch_pvc_numeric_by <- function() {
  f <- gamlss:::gamlss.pvc
  src <- deparse(f, width.cutoff = 500L)
  if (any(grepl("pvc_ctr", src, fixed = TRUE))) return(invisible(TRUE)) #already patched (file sourced again)
  old_gaic <-"GAIC <- sum(w * (y - fv)^2) + k * fit$edf"
  old_pred <- "pred <- drop(nx %*% fit$beta) * by.varExtra"
  i_gaic <- grep(old_gaic, src, fixed = TRUE)
  i_pred <- grep(old_pred, src, fixed = TRUE)
  i_by <- grep("by.var <- if (is.null(xeval))", src, fixed = TRUE)
  i_if <- grep("^\\s*if \\(!is.null\\(by.var\\)\\) \\{$", src)
  i_if <- i_if[i_if > i_by[1]][1]
  if (length(i_gaic) != 1 || length(i_pred) != 1 || length(i_by) != 1 || is.na(i_if))
    stop("patch_pvc_numeric_by: gamlss.pvc source not as expected (gamlss ", packageVersion("gamlss"), "); not fitting with an unpatched pvc")
  src[i_gaic] <- sub(old_gaic, "GAIC <- sum(w * (yvar - fv)^2) + k * fit$edf", src[i_gaic], fixed = TRUE)
  src[i_pred] <- sub(old_pred, "pred <- drop(nx %*% fit$beta) * (by.varExtra - pvc_ctr)", src[i_pred], fixed = TRUE)
  src <- append(src, "    pvc_ctr <- if (!is.null(by.var) && !is.factor(by.var)) mean(by.var) else 0; if (pvc_ctr != 0) by.var <- by.var - pvc_ctr",
                after = i_if - 1)
  g <- eval(parse(text = src))
  environment(g) <- environment(f)
  assignInNamespace("gamlss.pvc", g, "gamlss")
  invisible(TRUE)
}
patch_pvc_numeric_by()
#summarise age (in days) & sex of the df a model was fit on, for model_sums/*_summary.csv
#age back-transformed from logAge_days (log10(age_days)); falls back to age_days for models fit on linear age
df_demo_stats <- function(df) {
  if ("logAge_days" %in% names(df)) {
    age <- 10^df$logAge_days
  } else {
    age <- df$age_days
  }
  sex_male <- as.numeric(as.character(df$sexMale))
  data.frame(
    "age_days_min" = min(age),
    "age_days_mean" = mean(age),
    "age_days_max" = max(age),
    "age_days_sd" = sd(age),
    "pct_female" = 100 * mean(sex_male == 0)
  )
}
