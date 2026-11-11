#check_patch_pvc_gaic.R
#Validates patch_pvc_gaic() (code/gamlss_fit_funs.R): that the fix is needed on this gamlss and that it is right.
#
#The bug: for pvc(x, by = <numeric>), gamlss.pvc fits beta(x) to yvar = y / by with weights w * by^2, but its GAIC
#objective (fnGAIC) scores sum(w * (y - fv)^2), i.e. the raw working response against beta(x). With method = "GAIC"
#nlminb then runs lambda to the 1e7 cap (a linear sex x age term) whatever the data say.
#
#Checks
#  1 source:     the patch rewrites exactly one line (fnGAIC's), is idempotent, and refuses an unexpected source
#  2 necessary:  unpatched, pvc lambda hits the 1e7 cap both with and without a real (nonlinear) sex x age effect
#  3 correct:    on the exact inputs of the final backfitting call, patched gamlss.pvc's lambda minimizes an
#                independent implementation of the intended criterion, sum(w by^2 (yvar - fv)^2) + k edf;
#                the unpatched lambda (1e7) does not
#  4 sensible:   patched, a real nonlinear effect gives edf > 2, a better fit to the true mu and a lower BIC; no
#                effect gives a ~linear term (doesn't overfit)
#  5 unaffected: pvc with a factor by, and pvc without by, fit identically patched vs unpatched
#
#Needs only gamlss (no gamlssTools). Run from the repo root:
#  Rscript code/checks/check_patch_pvc_gaic.R
#FIT_FUNS overrides the path to gamlss_fit_funs.R; N_SIM the sample size.

suppressMessages(library(gamlss))

fit_funs <- Sys.getenv("FIT_FUNS", "code/gamlss_fit_funs.R")
n <- as.integer(Sys.getenv("N_SIM", 2000))
set.seed(20261108)

cat("gamlss", as.character(packageVersion("gamlss")), "\n")

#keep the unpatched smoother before anything patches it
pvc_orig <- gamlss:::gamlss.pvc
set_pvc <- function(f) assignInNamespace("gamlss.pvc", f, "gamlss")

#evaluate only patch_pvc_gaic()'s definition from gamlss_fit_funs.R (sourcing the whole file would patch on load
#and require gamlssTools)
exprs <- parse(fit_funs, keep.source = FALSE)
is_def <- vapply(exprs, function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                   identical(e[[2]], as.name("patch_pvc_gaic")), logical(1))
stopifnot("patch_pvc_gaic() definition not found in FIT_FUNS" = sum(is_def) == 1)
eval(exprs[[which(is_def)]])

results <- list()
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  results[[name]] <<- ok
  cat(sprintf("[%s] %s%s\n", if (ok) "PASS" else "FAIL", name, if (nzchar(detail)) paste0(" - ", detail) else ""))
}
src_of <- function(f) deparse(f, width.cutoff = 500L)

### 1 SOURCE ###
cat("\n## 1 source\n")
set_pvc(pvc_orig)
patch_pvc_gaic()
pvc_patched <- gamlss:::gamlss.pvc
s0 <- src_of(pvc_orig)
s1 <- src_of(pvc_patched)
changed <- which(s0 != s1)
check("patch changes exactly one line", length(s0) == length(s1) && length(changed) == 1,
      paste(length(changed), "line(s) differ"))
fn_start <- grep("fnGAIC <- function(lambda, k)", s0, fixed = TRUE)
fn_end <- grep("fnGAICforFactor <- function(lambda, k)", s0, fixed = TRUE)
check("changed line is in fnGAIC (numeric-by objective)",
      length(changed) == 1 && length(fn_start) == 1 && changed > fn_start && changed < fn_end)
check("changed line now scores yvar",
      length(changed) == 1 && grepl("sum(w * (yvar - fv)^2)", s1[changed], fixed = TRUE),
      if (length(changed) == 1) trimws(s1[changed]) else "")
check("patched function keeps the gamlss namespace as its environment",
      identical(environment(pvc_patched), environment(pvc_orig)))
patch_pvc_gaic() #second call (gamlss_fit_funs.R sourced twice in one session)
check("patch is idempotent", identical(src_of(gamlss:::gamlss.pvc), s1))

#a gamlss.pvc that has neither the old nor twice the new line (e.g. another gamlss version) must stop the fit
bad <- pvc_orig
body(bad) <- str2lang(paste(sub("GAIC <- sum(w * (y - fv)^2) + k * fit$edf", "GAIC <- sum(w * (y - fv)^2)",
                                s0[-1], fixed = TRUE), collapse = "\n"))
environment(bad) <- environment(pvc_orig)
set_pvc(bad)
err <- tryCatch({patch_pvc_gaic(); NULL}, error = function(e) conditionMessage(e))
check("patch refuses an unexpected gamlss.pvc source", !is.null(err), if (!is.null(err)) err else "no error")
set_pvc(pvc_orig)

### SIMULATED DATA ###
#log-age-like x on [0, 1], unequal sex sizes so pvc's centring of by is not a no-op
x <- runif(n)
sexEffect <- ifelse(runif(n) < 0.55, 1, -1)
sex <- factor(ifelse(sexEffect == 1, "M", "F"))
f_true <- function(x) 2 + 1.5 * plogis((x - 0.3) / 0.08)
g_true <- function(x) 0.25 * sin(2 * pi * x) #half the M - F gap: clearly nonlinear
sigma <- 0.3
mu_sig <- f_true(x) + sexEffect * g_true(x)
mu_null <- f_true(x)
dat <- data.frame(x, sexEffect, sex,
                  y_sig = mu_sig + rnorm(n, 0, sigma),
                  y_null = mu_null + rnorm(n, 0, sigma))
k <- log(n)

#the pipeline's mu term: pb(age) + pvc(age, by = sexEffect) + sexEffect, lambdas by GAIC with k = log(n)
fit_pvc <- function(resp) {
  f <- as.formula(paste(resp, "~ pb(x, method = 'GAIC', k = k) + pvc(x, by = sexEffect, method = 'GAIC', k = k) + sexEffect"))
  gamlss(f, sigma.formula = ~1, family = NO, data = dat, trace = FALSE,
         control = gamlss.control(n.cyc = 200))
}
pvc_term <- function(m) {
  s <- m$mu.coefSmo[[which(vapply(m$mu.coefSmo, inherits, logical(1), "pvc"))]]
  c(lambda = s$lambda, edf = s$edf)
}
rmse <- function(a, b) sqrt(mean((a - b)^2))
at_cap <- function(l) l >= 1e7 * (1 - 1e-8)

#record the inputs of the last (final backfitting) call to gamlss.pvc
rec <- new.env()
recording <- function(fun) {
  function(x, y, w, xeval = NULL, ...) {
    if (is.null(xeval)) { rec$x <- x; rec$y <- y; rec$w <- w }
    fun(x, y, w, xeval, ...)
  }
}

set_pvc(recording(pvc_orig))
m_orig_sig <- fit_pvc("y_sig")
in_orig <- as.list(rec)
m_orig_null <- fit_pvc("y_null")
set_pvc(recording(pvc_patched))
m_pat_sig <- fit_pvc("y_sig")
in_pat <- as.list(rec)
m_pat_null <- fit_pvc("y_null")
set_pvc(pvc_orig)
stopifnot("gamlss.pvc wrapper was never called" = !is.null(in_pat$x))

t_os <- pvc_term(m_orig_sig); t_on <- pvc_term(m_orig_null)
t_ps <- pvc_term(m_pat_sig); t_pn <- pvc_term(m_pat_null)
cat("\npvc term (lambda, edf):\n")
print(rbind(unpatched_signal = t_os, unpatched_null = t_on, patched_signal = t_ps, patched_null = t_pn))

### 2 NECESSARY ###
cat("\n## 2 necessary\n")
check("unpatched: lambda at 1e7 cap with a real nonlinear effect", at_cap(t_os["lambda"]),
      sprintf("lambda = %.3g, edf = %.2f", t_os["lambda"], t_os["edf"]))
check("unpatched: lambda at 1e7 cap with no effect", at_cap(t_on["lambda"]),
      sprintf("lambda = %.3g", t_on["lambda"]))

### 3 CORRECT ###
cat("\n## 3 correct\n")
#independent implementation of the intended criterion on one gamlss.pvc call's inputs
reference <- function(inp) {
  X <- as.matrix(attr(inp$x, "X")); D <- as.matrix(attr(inp$x, "D")); by <- attr(inp$x, "by")
  kk <- attr(inp$x, "control")$k
  W <- inp$w * by^2
  yvar <- inp$y / by
  XtWX <- crossprod(X, W * X); XtWy <- crossprod(X, W * yvar); DtD <- crossprod(D)
  crit <- function(lambda) {
    A <- XtWX + lambda * DtD
    beta <- solve(A, XtWy)
    edf <- sum(diag(solve(A, XtWX)))
    sum(W * (yvar - X %*% beta)^2) + kk * edf
  }
  #coarse grid on log10(lambda), then refine around the best point
  lg <- seq(-7, 7, by = 0.1)
  cg <- vapply(10^lg, crit, numeric(1))
  i <- which.min(cg)
  opt <- optimize(function(l) crit(10^l), c(lg[max(1, i - 1)], lg[min(length(lg), i + 1)]), tol = 1e-8)
  list(crit = crit, lambda = 10^opt$minimum, min = opt$objective, by_centred = abs(mean(by)) < 1e-10,
       by_raw = all(abs(by - (sexEffect - mean(sexEffect))) < 1e-10))
}
#rerun one gamlss.pvc call on recorded inputs from a fixed lambda start (default pvc start)
rerun <- function(fun, inp, start = 1) {
  assign(as.character(attr(inp$x, "NameForLambda")), start, envir = attr(inp$x, "gamlss.env"))
  fun(inp$x, inp$y, inp$w)$lambda
}

ref <- reference(in_pat)
check("recorded by is sexEffect centred on its mean", ref$by_raw && ref$by_centred)
#starts: pvc.control()'s default (10) and the lambdas backfitting passes on between calls in a real fit
for (st in c(pvc.control()$start, 1, 1e3, unname(t_ps["lambda"]))) {
  l_pat <- rerun(pvc_patched, in_pat, st)
  gap <- ref$crit(l_pat) - ref$min
  check(sprintf("patched lambda minimizes the intended GAIC (start %.4g)", st),
        gap < 1e-6 * abs(ref$min) + 1e-6,
        sprintf("patched %.4g vs reference %.4g; criterion gap %.2e", l_pat, ref$lambda, gap))
}
#not a patch check: gamlss.pvc's nlminb works on raw lambda, so from a start far above the optimum (where the
#criterion is nearly flat) it can stop early. Reported for information only
l_far <- rerun(pvc_patched, in_pat, 1e6)
cat(sprintf("[info] patched, start 1e6: lambda %.4g (criterion gap %.2e) - nlminb start dependence, not the patch\n",
            l_far, ref$crit(l_far) - ref$min))
check("fit's own pvc lambda matches the reference", abs(log10(t_ps["lambda"] / ref$lambda)) < 0.05,
      sprintf("%.4g vs %.4g", t_ps["lambda"], ref$lambda))
l_unp <- rerun(pvc_orig, in_pat, 1)
check("unpatched lambda on the same inputs is the cap and scores worse on the intended GAIC",
      at_cap(l_unp) && ref$crit(l_unp) - ref$min > 1,
      sprintf("lambda %.3g; criterion %.1f vs minimum %.1f", l_unp, ref$crit(l_unp), ref$min))
ref_o <- reference(in_orig)
check("on the unpatched fit's own final inputs, the intended GAIC also prefers lambda < cap",
      ref_o$lambda < 1e6, sprintf("reference lambda %.4g", ref_o$lambda))

### 4 SENSIBLE ###
cat("\n## 4 sensible\n")
r_o <- rmse(fitted(m_orig_sig), mu_sig); r_p <- rmse(fitted(m_pat_sig), mu_sig)
check("patched: nonlinear effect gives edf > 2", t_ps["edf"] > 2.5, sprintf("edf %.2f", t_ps["edf"]))
check("patched: closer to the true mu than unpatched", r_p < r_o, sprintf("RMSE %.4f vs %.4f", r_p, r_o))
b_o <- GAIC(m_orig_sig, k = k); b_p <- GAIC(m_pat_sig, k = k)
check("patched: lower BIC than unpatched", b_p < b_o, sprintf("BIC %.2f vs %.2f", b_p, b_o))
check("patched: no effect gives a ~linear term (edf < 3)", t_pn["edf"] < 3,
      sprintf("lambda %.3g, edf %.2f", t_pn["lambda"], t_pn["edf"]))

### 5 UNAFFECTED ###
cat("\n## 5 unaffected\n")
fit_both <- function(f) {
  fits <- lapply(list(pvc_orig, pvc_patched), function(fun) {
    set_pvc(fun)
    m <- gamlss(f, sigma.formula = ~1, family = NO, data = dat, trace = FALSE,
                control = gamlss.control(n.cyc = 200))
    c(deviance = deviance(m), m$mu.lambda)
  })
  set_pvc(pvc_orig)
  fits
}
fac <- fit_both(y_sig ~ sex + pvc(x, by = sex, method = "GAIC", k = k))
check("pvc(x, by = factor): identical patched vs unpatched", isTRUE(all.equal(fac[[1]], fac[[2]], tolerance = 1e-10)))
noby <- fit_both(y_sig ~ pvc(x, method = "GAIC", k = k))
check("pvc(x) without by: identical patched vs unpatched", isTRUE(all.equal(noby[[1]], noby[[2]], tolerance = 1e-10)))

### SUMMARY ###
n_fail <- sum(!unlist(results))
cat(sprintf("\n%d/%d checks passed\n", length(results) - n_fail, length(results)))
if (n_fail > 0) stop(n_fail, " check(s) failed")
