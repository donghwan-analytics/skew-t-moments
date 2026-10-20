############################################################
## [STANDARD full skew-t, alpha != 0]
## Compare:
##   (A) my          : handwritten standard-form skew-t raw moments
##   (B) exact_numint: numerical integration of the actual skew-t density
## Optional:
##   (C) mc_check    : Monte Carlo check using sn::rmst
##
## Target distribution:
##   T ~ Skew-t_2(0, I, alpha, nu)
##
## Handwritten formula implemented exactly:
##   T = Lambda^{-1/2} Z,
##   Lambda ~ Gamma(nu/2, nu/2), independent of Z,
##   Z ~ SN_2(0, I, alpha)
##
##   Z = delta |U| + Y,
##   U ~ N(0,1), Y ~ N_2(0, Sigma), U independent of Y,
##   delta = alpha / sqrt(1 + alpha^T alpha),
##   Sigma = I - delta delta^T.
##
##   For k = (k1,k2), |k| = k1 + k2 < nu,
##
##   E[T^k] = E[Lambda^{-|k|/2}] E[Z^k]
##
##   E[Z^k] = sum_{0 <= r <= k}
##              [ prod_i choose(k_i, r_i) delta_i^{k_i-r_i} ]
##              E|U|^{|k|-|r|} E[Y^r].
############################################################

## ------------------------------
## 0) package helper
## ------------------------------
load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

load_or_install("sn")
load_or_install("cubature")

## ------------------------------
## 1) utilities
## ------------------------------
E_Lambda_neg_half <- function(m, nu) {
  if (m >= nu) return(NA_real_)
  (nu / 2)^(m / 2) * gamma((nu - m) / 2) / gamma(nu / 2)
}

E_absU_power <- function(p) {
  if (p <= -1) stop("Need p > -1 for E|U|^p.")
  2^(p / 2) * gamma((p + 1) / 2) / sqrt(pi)
}

isserlis_moment <- function(Sigma, idx) {
  m <- length(idx)
  if (m == 0) return(1)
  if (m %% 2 == 1) return(0)

  i1 <- idx[1]
  total <- 0
  for (j in 2:m) {
    cov_ij <- Sigma[i1, idx[j]]
    rest <- idx[-c(1, j)]
    total <- total + cov_ij * isserlis_moment(Sigma, rest)
  }
  total
}

E_Y_moment <- function(r1, r2, Sigma) {
  idx <- c(rep(1, r1), rep(2, r2))
  isserlis_moment(Sigma, idx)
}

## ------------------------------
## 2) my handwritten standard-form formula
## ------------------------------
my_standard_skewn_raw_moment_biv <- function(k1, k2, alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 2) stop("alpha must be length 2.")

  delta <- alpha / sqrt(1 + sum(alpha^2))
  Sigma <- diag(2) - tcrossprod(delta)

  out <- 0
  for (r1 in 0:k1) {
    for (r2 in 0:k2) {
      p <- (k1 + k2) - (r1 + r2)

      coef_part <- choose(k1, r1) * choose(k2, r2) *
        delta[1]^(k1 - r1) * delta[2]^(k2 - r2)

      out <- out + coef_part * E_absU_power(p) * E_Y_moment(r1, r2, Sigma)
    }
  }
  out
}

my_standard_skewt_raw_moment_biv <- function(k1, k2, alpha, nu) {
  m <- k1 + k2
  if (m >= nu) return(NA_real_)
  E_Lambda_neg_half(m, nu) * my_standard_skewn_raw_moment_biv(k1, k2, alpha)
}

## ------------------------------
## 3) actual skew-t moment by numerical integration
##    exact_numint := integral against sn::dmst density
## ------------------------------
exact_skewt_raw_moment_numint_biv <- function(k1, k2, alpha, nu,
                                              tol = 1e-7,
                                              maxEval = 200000,
                                              eps = 1e-12) {
  m <- k1 + k2
  if (m >= nu) return(NA_real_)

  xi <- c(0, 0)
  Omega <- diag(2)

  integrand_u <- function(u) {
    x1 <- tan(pi * (u[1] - 0.5))
    x2 <- tan(pi * (u[2] - 0.5))
    x  <- c(x1, x2)

    jac <- pi^2 * (1 + x1^2) * (1 + x2^2)

    dens <- as.numeric(sn::dmst(x = x, xi = xi, Omega = Omega,
                                alpha = alpha, nu = nu))

    (x1^k1) * (x2^k2) * dens * jac
  }

  out <- cubature::hcubature(
    f          = integrand_u,
    lowerLimit = c(eps, eps),
    upperLimit = c(1 - eps, 1 - eps),
    tol        = tol,
    maxEval    = maxEval,
    vectorInterface = FALSE
  )

  as.numeric(out$integral)
}

## ------------------------------
## 4) optional Monte Carlo check
## ------------------------------
mc_skewt_raw_moment_biv <- function(k1, k2, alpha, nu,
                                    n = 300000,
                                    seed = 123) {
  set.seed(seed)
  X <- sn::rmst(n = n, xi = c(0, 0), Omega = diag(2), alpha = alpha, nu = nu)
  vals <- X[, 1]^k1 * X[, 2]^k2
  c(mean = mean(vals), se = sd(vals) / sqrt(n), sd = sd(vals))
}

## ------------------------------
## 5) comparison runner
## ------------------------------
compare_standard_skewt <- function(alpha, nu,
                                   moment_list,
                                   do_mc = FALSE,
                                   mc_n = 300000,
                                   tol = 1e-7,
                                   maxEval = 200000) {
  out <- data.frame()

  for (kk in moment_list) {
    k1 <- kk[1]
    k2 <- kk[2]
    ord <- k1 + k2
    if (ord >= nu) next

    my_val <- my_standard_skewt_raw_moment_biv(k1, k2, alpha, nu)
    ex_val <- exact_skewt_raw_moment_numint_biv(k1, k2, alpha, nu,
                                                tol = tol, maxEval = maxEval)

    one <- data.frame(
      alpha1 = alpha[1],
      alpha2 = alpha[2],
      nu     = nu,
      k1     = k1,
      k2     = k2,
      order  = ord,
      my     = my_val,
      exact_numint = ex_val,
      diff    = my_val - ex_val,
      abs_diff = abs(my_val - ex_val)
    )

    if (do_mc) {
      mc <- mc_skewt_raw_moment_biv(k1, k2, alpha, nu, n = mc_n,
                                    seed = 1000 + 100 * k1 + k2)
      one$mc_mean <- unname(mc["mean"])
      one$mc_se   <- unname(mc["se"])
      one$my_minus_mc <- my_val - one$mc_mean
      one$exact_minus_mc <- ex_val - one$mc_mean
    }

    out <- rbind(out, one)
  }

  out[order(out$order, out$k1, out$k2), ]
}

## ------------------------------
## 6) moment list and parameter grid
## ------------------------------
moment_list <- list(
  c(1, 0), c(0, 1), c(1, 1),
  c(2, 0), c(0, 2), c(2, 1), c(1, 2),
  c(3, 0), c(0, 3),
  c(2, 2), c(4, 0), c(0, 4)
)

param_grid <- list(
  list(name = "std_case_1", alpha = c( 0.8, -1.2), nu =  8),
  list(name = "std_case_2", alpha = c( 1.5,  0.7), nu = 10),
  list(name = "std_case_3", alpha = c(-2.0,  1.0), nu = 12),
  list(name = "std_case_4", alpha = c( 0.5,  2.5), nu = 15)
)

## ------------------------------
## 7) run
## ------------------------------
res <- do.call(
  rbind,
  lapply(param_grid, function(pg) {
    ans <- compare_standard_skewt(
      alpha   = pg$alpha,
      nu      = pg$nu,
      moment_list = moment_list,
      do_mc   = FALSE,   # set TRUE if you want MC check too
      mc_n    = 300000,
      tol     = 1e-7,
      maxEval = 200000
    )
    ans$case <- pg$name
    ans
  })
)

res <- res[, c("case", "alpha1", "alpha2", "nu", "k1", "k2", "order",
               "my", "exact_numint", "diff", "abs_diff")]

cat("\n================ FULL RESULT TABLE ================\n")
print(res, row.names = FALSE)

cat("\n================ OVERALL MAX ABS DIFF ================\n")
cat("max_abs_diff = ", max(res$abs_diff, na.rm = TRUE), "\n", sep = "")

cat("\n================ ONLY ORDER 4 ================\n")
print(subset(res, order == 4), row.names = FALSE)

cat("\n================ SUMMARY BY CASE ================\n")
summary_by_case <- do.call(
  rbind,
  lapply(split(res, res$case), function(df) {
    data.frame(
      case = df$case[1],
      max_abs_diff = max(df$abs_diff, na.rm = TRUE)
    )
  })
)
print(summary_by_case, row.names = FALSE)
