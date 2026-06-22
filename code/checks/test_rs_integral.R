# Tests for the updated rs_integral function.
# Run from the project root: Rscript code/test_rs_integral.R
#
# Verifies that rs_integral(logAge, y, g = function(x) 10^x) returns
# Sum y_i * Delta(age_i), i.e. a left Riemann-Stieltjes approximation
# of the integral of y over raw age, given y sampled at log-spaced age points.

rs_integral <- function(x, y, g = function(x) 10^(x)) {
  stopifnot(length(x) == length(y), length(x) >= 2)
  ord <- order(x); x <- x[ord]; y <- y[ord]
  gx <- g(x)
  if (any(diff(gx) < 0)) stop("g(x) must be non-decreasing for a Riemann-Stieltjes integral")
  sum(y[-length(y)] * diff(gx))
}

make_grid <- function(a, b, n) {
  logAge <- seq(log10(a), log10(b), length.out = n)
  list(logAge = logAge, age = 10^logAge)
}

run_test <- function(label, fy, a, b, expected, n = 1000, tol = NULL, exact = FALSE) {
  grid <- make_grid(a, b, n)
  y <- fy(grid$age)
  est <- rs_integral(grid$logAge, y)
  err <- est - expected
  rel <- err / expected
  if (exact) {
    stopifnot(abs(err) < 1e-8)
    cat(sprintf("[EXACT  ] %-40s est=%.6g expected=%.6g abs_err=%.2e\n",
                label, est, expected, abs(err)))
  } else {
    cat(sprintf("[APPROX ] %-40s n=%-5d est=%.6g expected=%.6g rel_err=%+.3e\n",
                label, n, est, expected, rel))
    if (!is.null(tol)) stopifnot(abs(rel) < tol)
  }
  invisible(est)
}

# 1. Constant y -- exact regardless of grid spacing
run_test("constant y=1, age in [1,100]",
         function(a) rep(1, length(a)),
         1, 100, expected = 99, exact = TRUE)

# 2. y = age -- converges as n grows
run_test("y=age, [1,100], n=100",   function(a) a, 1, 100, expected = (100^2 - 1^2)/2, n = 100)
run_test("y=age, [1,100], n=1000",  function(a) a, 1, 100, expected = (100^2 - 1^2)/2, n = 1000)
run_test("y=age, [1,100], n=10000", function(a) a, 1, 100, expected = (100^2 - 1^2)/2, n = 10000)

# 3. y = age^2
run_test("y=age^2, [1,100], n=10000",
         function(a) a^2, 1, 100,
         expected = (100^3 - 1^3)/3, n = 10000, tol = 1e-2)

# 4. y = 1/age
run_test("y=1/age, [1,100], n=10000",
         function(a) 1/a, 1, 100,
         expected = log(100), n = 10000, tol = 1e-2)

# 5. Constant y on realistic LBCC age range -- exact
run_test("constant y=-2, age in days [280, 36805]",
         function(a) rep(-2, length(a)),
         280, 280 + 100*365.25,
         expected = -2 * (100*365.25), exact = TRUE)

# 6a. Cross-check vs pracma::trapz on a smooth, non-oscillating integrand.
# rs_integral is a left Riemann-Stieltjes sum, so it lags trapezoidal by O(Delta age);
# on a smooth y with n=10000 the relative gap should be small (<1%).
if (requireNamespace("pracma", quietly = TRUE)) {
  grid <- make_grid(1, 100, 10000)
  y <- 1 / (1 + grid$age)  # smooth, monotone, non-oscillating
  rsi <- rs_integral(grid$logAge, y)
  trp <- pracma::trapz(grid$age, y)
  cat(sprintf("[CROSS  ] %-40s rs_integral=%.6g trapz=%.6g rel_diff=%+.3e\n",
              "y=1/(1+age), [1,100], n=10000",
              rsi, trp, (rsi - trp) / trp))
  stopifnot(abs(rsi - trp) / abs(trp) < 1e-2)

  # 6b. With an oscillating integrand the absolute gap at any fixed n can be
  # large relative to a small true integral, but |rs - trapz| must shrink as n
  # grows -- that's the property that pins down rs_integral as a left Riemann
  # sum approximating the same integral.
  ns <- c(1e3, 1e4, 1e5)
  diffs <- sapply(ns, function(n) {
    g <- make_grid(1, 100, n)
    abs(rs_integral(g$logAge, sin(g$age)) - pracma::trapz(g$age, sin(g$age)))
  })
  cat(sprintf("[CONVRG ] sin(age) |rs - trapz| at n=%s: %s\n",
              paste(ns, collapse=", "),
              paste(sprintf("%.3e", diffs), collapse=", ")))
  stopifnot(diffs[2] < diffs[1], diffs[3] < diffs[2])
} else {
  cat("[SKIP   ] pracma not installed -- skipping trapz cross-checks\n")
}

# Sign/orientation: flipping y flips the integral
grid <- make_grid(1, 100, 1000); y <- grid$age
stopifnot(all.equal(rs_integral(grid$logAge, -y), -rs_integral(grid$logAge, y)))
cat("[OK     ] sign flip\n")

# Monotonicity guard rejects a decreasing integrator
err <- tryCatch(rs_integral(grid$logAge, y, g = function(x) -10^x),
                error = function(e) e)
stopifnot(inherits(err, "error"))
cat("[OK     ] monotonicity guard rejects decreasing g\n")

# Discrete identity: rs_integral with g(x)=10^x must equal a direct left
# Riemann sum of y over the age grid (since gx == age).
manual_left <- sum(head(y, -1) * diff(grid$age))
stopifnot(all.equal(rs_integral(grid$logAge, y), manual_left))
cat("[OK     ] matches direct left Riemann sum on age grid\n")

# Change-of-variables consistency (asymptotic, NOT exact at finite n):
# Integral y d(age) = Integral y * age * ln(10) d(logAge), but the two
# left Riemann sums of those integrands differ by O(Delta logAge^2) per term
# (Taylor: age_{i+1} - age_i = age_i * ln(10) * d(logAge) * (1 + O(d(logAge)))).
# So we verify the gap shrinks as n grows.
chvar_rel_diff <- function(n) {
  g <- make_grid(1, 100, n)
  yy <- g$age
  rsi <- rs_integral(g$logAge, yy)
  manual <- sum(head(yy, -1) * head(g$age, -1) * log(10) * diff(g$logAge))
  abs(rsi - manual) / abs(rsi)
}
ns_chvar <- c(1e3, 1e4, 1e5)
diffs_chvar <- sapply(ns_chvar, chvar_rel_diff)
cat(sprintf("[CONVRG ] change-of-variables rel gap at n=%s: %s\n",
            paste(formatC(ns_chvar, format="d"), collapse=", "),
            paste(sprintf("%.3e", diffs_chvar), collapse=", ")))
stopifnot(diffs_chvar[2] < diffs_chvar[1], diffs_chvar[3] < diffs_chvar[2])

cat("\nAll tests passed.\n")
