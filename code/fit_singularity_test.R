#Fit GAMLSS models to select from on CV samples

#LOAD PACKAGES
library(data.table)
library(dplyr)
library(gamlss)
library(gamlssTools)

base <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/"
source(paste0(base, "code/gamlss_fit_funs.R"))
options(warn = 1)

# ################
# #testing new safe_gamlss version
# safe_gamlss <- function (...) 
# {
#   print("new safe_gamlss")
#   warn_list <- NULL
#   mod <- withCallingHandlers({
#     gamlss(...)
#   }, warning = function(w) {
#     warn_list <<- c(warn_list, w$message)
#     if (grepl("Error", w$message, ignore.case = TRUE) || 
#         grepl("converge", w$message, ignore.case = TRUE)) {
#       stop(simpleError(w$message))
#     }
#   }, error = function(e) {
#     stop(e)
#   })
#   null_mu <- is.null(coef(mod, what = "mu"))
#   null_sigma <- is.null(coef(mod, what = "sigma"))
#   if (null_mu && null_sigma) {
#     stop("Model fit failed: coefficients are NULL")
#   }
#   if (mod$converged == FALSE) {
#     stop("Model did not converge:", warn_msg)
#   }
#   
#   #check warnings again
#   warn_combined <- paste(warn_list, collapse = " | ")
#   if (grepl("Error", warn_combined, ignore.case = TRUE) || 
#       grepl("converge", warn_combined, ignore.case = TRUE)) {
#     stop(paste("Model failed due to warning:", warn_combined))
#   }
#   
#   return(mod)
# }
# #####################

warning("this is a warning")

df_path <- "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/data/cv_sample_B_dfs/SUBC.Right.Lateral.Ventricle_totalFALSE_logPhenoFALSE_logAgeTRUE.csv"                                                       

#GET ARGS
args <- commandArgs(trailingOnly = TRUE)
gamlss_form <- args[1]
print(gamlss_form)
df <- fread(df_path, stringsAsFactors = TRUE, na.strings = "") #path to csv

#FIT BASIC MODEL
print(gamlss_form)
model <- eval(parse(text = gamlss_form))

saveRDS(model, "/mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/model.rds")

print("SUCCESS")
