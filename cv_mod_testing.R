#troubleshooting
# gamlss_3lambda_iris <- function(pheno, lambda=NULL, 
#                            fam="GG",
#                            nu_form="1",
#                            start.from=NULL){
#   
# 
#   #define formulas for each moment
#   mu_form <- paste(
#     "gamlss(formula =", pheno, "~ pb(Petal.Length) + Species")
#   
#   
#   sig_form <- paste(
#     "sigma.formula = ~ Petal.Length"
#   )
#   
#   nu_form <- paste("nu.formula = ~", nu_form)
#   
#   if (is.null(start.from)) {
#     control <- paste("control = gamlss.control(n.cyc = 200), family =", fam, ", data= iris, trace = FALSE)")
#   } else if (is.gamlss(start.from)) {
#     control <- paste0("start.from = ", start.from,
#                       ", control = gamlss.control(n.cyc = 200), family =", fam,
#                       ", data= iris, trace = FALSE)")
#   } else {
#     stop("start.from arg must be gamlss model")
#   }
# 
#   #try methods
#   
#   result <- tryCatch({
#     gamlss_RSformula <-paste(mu_form, sig_form, nu_form, control, sep=",")
#     print(gamlss_RSformula)
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
#       gamlss_CGformula <-paste(mu_form, sig_form, nu_form, "method=CG()", control, sep=",")
#       eval(parse(text = gamlss_CGformula))
#       
#       #if CG also fails, return NULL
#     }, error = function(e2) {
#       message(e2$message, ", returning NULL")
#       return(NULL)
#     })
#   } , finally = {
#     message("done")
#   } )
#   return(result)
# }

gamlss_3lambda_iris <- function(pheno, lambda=NULL, 
                                fam="GG",
                                nu_form="1",
                                start.from=NULL) {
  
  mu_formula <- as.formula(paste(pheno, "~ pb(Petal.Length) + Species"))
  sigma_formula <- ~ Petal.Length
  nu_formula <- as.formula(paste("~", nu_form))
  
  args <- list(
    formula = mu_formula,
    sigma.formula = sigma_formula,
    nu.formula = nu_formula,
    family = fam,
    data = iris,
    control = gamlss.control(n.cyc = 200),
    trace = FALSE
  )
  
  if (!is.null(start.from)) {
    if (!is.gamlss(start.from)) stop("start.from arg must be gamlss model")
    args$start.from <- start.from
  }
  
  result <- tryCatch({
    do.call(gamlss, args)
  }, warning = function(w) {
    message("warning: ", conditionMessage(w))
    do.call(gamlss, args)
  }, error = function(e) {
    message("error: ", conditionMessage(e), ", trying method = CG()")
    args$method <- CG()
    tryCatch({
      do.call(gamlss, args)
    }, error = function(e2) {
      message("error: ", conditionMessage(e2), ", returning NULL")
      return(NULL)
    })
  }, finally = {
    message("done")
  })
  
  return(result)
}


nu_list <- list(int = "1", 
                species = "Species", 
                plength = "Petal.Length")

mod_list <- c()

for (nu in nu_list){
  nu_name <- names(nu_list)[nu_list==nu]
  print(nu_name)
  
  model <- gamlss_3lambda_iris("Sepal.Width", nu_form=nu, start.from=mod_list[[1]])
  
  #if model isn't fit, skip to next loop
  if (is.null(model)) {
    message("model fitting failed")
    next
  }
  
  mod_list[[nu_name]] <- model
}
  