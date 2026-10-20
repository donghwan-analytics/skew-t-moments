############################################################
## [GENERAL full skew-t, alpha != 0]
## Compare:
##   (A) my           : handwritten general-form skew-t raw moments
##   (B) exact_numint : numerical integration of the actual skew-t density
## Optional:
##   (C) mc_check     : Monte Carlo check using sn::rmst
##
## Target distribution:
##   Y ~ Skew-t_2(xi, Omega, alpha, nu)
##
## Handwritten formula implemented exactly:
##   omega = diag(sqrt(Omega_11), sqrt(Omega_22))
##   barOmega = omega^{-1} Omega omega^{-1}
##   choose symmetric square root B = barOmega^{1/2}
##   A = omega B,  so that  A A^T = Omega
##   beta = B alpha
##
##   T = A^{-1}(Y - xi) ~ ST_2(0, I, beta, nu)
##   Y = xi + A T
##
##   For k = (k1,k2),
##   E[(Y-xi)^k]
##     = sum_{N in U(k)}
##         [ prod_i k_i! / prod_j n_ij! ]
##         [ prod_i prod_j a_ij^{n_ij} ]
##         M_std_{c(N)}(beta, nu)
##
##   where c(N) is the column-sum multi-index,
##   M_std_r(beta, nu) = E[T^r].
##
##   Then raw moment:
##   E[Y^k] = sum_{s <= k} choose(k,s) xi^{k-s} E[(Y-xi)^s].
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

sym_sqrt <- function(M) {
  ee <- eigen(M, symmetric = TRUE)
  if (any(ee$values <= 0)) stop("Matrix must be positive definite.")
  ee$vectors %*% diag(sqrt(ee$values)) %*% t(ee$vectors)
}

## ------------------------------
## 2) my handwritten standard-form formula
##    M_std_r(beta, nu) = E[T^r]
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
## 3) enumerate N in U(k) for 2x2 case
##    N = [ [n11, n12],
##          [n21, n22] ]
##    with row sums (k1, k2)
## ------------------------------
enumerate_N_biv <- function(k1, k2) {
  out <- list()
  idx <- 1L
  for (n11 in 0:k1) {
    n12 <- k1 - n11
    for (n21 in 0:k2) {
      n22 <- k2 - n21
      out[[idx]] <- matrix(c(n11, n12,
                             n21, n22), nrow = 2, byrow = TRUE)
      idx <- idx + 1L
    }
  }
  out
}

## ------------------------------
## 4) centered moment from handwritten affine formula
## ------------------------------
my_general_skewt_centered_moment_biv <- function(k1, k2, Omega, alpha, nu) {
  m <- k1 + k2
  if (m >= nu) return(NA_real_)

  Omega <- as.matrix(Omega)
  alpha <- as.numeric(alpha)

  if (!all(dim(Omega) == c(2, 2))) stop("Omega must be 2x2.")
  if (length(alpha) != 2) stop("alpha must be length 2.")

  omega <- diag(sqrt(diag(Omega)))
  barOmega <- solve(omega) %*% Omega %*% solve(omega)
  B <- sym_sqrt(barOmega)
  A <- omega %*% B
  beta <- as.numeric(B %*% alpha)

  N_list <- enumerate_N_biv(k1, k2)
  out <- 0

  for (N in N_list) {
    n11 <- N[1, 1]; n12 <- N[1, 2]
    n21 <- N[2, 1]; n22 <- N[2, 2]

    coef_multinom <- factorial(k1) / (factorial(n11) * factorial(n12)) *
      factorial(k2) / (factorial(n21) * factorial(n22))

    coef_A <- (A[1,1]^n11) * (A[1,2]^n12) *
      (A[2,1]^n21) * (A[2,2]^n22)

    c1 <- n11 + n21
    c2 <- n12 + n22

    out <- out + coef_multinom * coef_A *
      my_standard_skewt_raw_moment_biv(c1, c2, beta, nu)
  }

  out
}

## ------------------------------
## 5) raw moment from handwritten binomial expansion
## ------------------------------
my_general_skewt_raw_moment_biv <- function(k1, k2, xi, Omega, alpha, nu) {
  m <- k1 + k2
  if (m >= nu) return(NA_real_)

  xi <- as.numeric(xi)
  out <- 0

  for (s1 in 0:k1) {
    for (s2 in 0:k2) {
      coef_part <- choose(k1, s1) * choose(k2, s2) *
        xi[1]^(k1 - s1) * xi[2]^(k2 - s2)

      out <- out + coef_part *
        my_general_skewt_centered_moment_biv(s1, s2, Omega, alpha, nu)
    }
  }

  out
}

## ------------------------------
## 6) actual skew-t moment by numerical integration
##    exact_numint := integral against sn::dmst density
## ------------------------------
exact_skewt_raw_moment_numint_biv <- function(k1, k2, xi, Omega, alpha, nu,
                                              tol = 1e-7,
                                              maxEval = 200000,
                                              eps = 1e-12) {
  m <- k1 + k2
  if (m >= nu) return(NA_real_)

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
## 7) optional Monte Carlo check
## ------------------------------
mc_skewt_raw_moment_biv <- function(k1, k2, xi, Omega, alpha, nu,
                                    n = 300000,
                                    seed = 123) {
  set.seed(seed)
  X <- sn::rmst(n = n, xi = xi, Omega = Omega, alpha = alpha, nu = nu)
  vals <- X[, 1]^k1 * X[, 2]^k2
  c(mean = mean(vals), se = sd(vals) / sqrt(n), sd = sd(vals))
}

## ------------------------------
## 8) comparison runner
## ------------------------------
compare_general_skewt <- function(xi, Omega, alpha, nu,
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

    my_val <- my_general_skewt_raw_moment_biv(k1, k2, xi, Omega, alpha, nu)
    ex_val <- exact_skewt_raw_moment_numint_biv(k1, k2, xi, Omega, alpha, nu,
                                                tol = tol, maxEval = maxEval)

    one <- data.frame(
      xi1     = xi[1],
      xi2     = xi[2],
      omega11 = Omega[1,1],
      omega12 = Omega[1,2],
      omega22 = Omega[2,2],
      alpha1  = alpha[1],
      alpha2  = alpha[2],
      nu      = nu,
      k1      = k1,
      k2      = k2,
      order   = ord,
      my      = my_val,
      exact_numint = ex_val,
      diff     = my_val - ex_val,
      abs_diff = abs(my_val - ex_val)
    )

    if (do_mc) {
      mc <- mc_skewt_raw_moment_biv(k1, k2, xi, Omega, alpha, nu,
                                    n = mc_n,
                                    seed = 2000 + 100 * k1 + k2)
      one$mc_mean <- unname(mc["mean"])
      one$mc_se   <- unname(mc["se"])
      one$my_minus_mc    <- my_val - one$mc_mean
      one$exact_minus_mc <- ex_val - one$mc_mean
    }

    out <- rbind(out, one)
  }

  out[order(out$order, out$k1, out$k2), ]
}

## ------------------------------
## 9) moment list and parameter grid
## ------------------------------
moment_list <- list(
  c(1, 0), c(0, 1), c(1, 1),
  c(2, 0), c(0, 2), c(2, 1), c(1, 2),
  c(3, 0), c(0, 3),
  c(2, 2), c(4, 0), c(0, 4)
)

param_grid <- list(
  list(
    name  = "gen_case_1",
    nu    = 10,
    xi    = c( 1.0, -0.5),
    Omega = matrix(c(2.0, 0.5,
                     0.5, 1.3), 2, 2, byrow = TRUE),
    alpha = c( 1.2, -0.8)
  ),
  list(
    name  = "gen_case_2",
    nu    = 12,
    xi    = c(-1.5,  2.0),
    Omega = matrix(c(1.4, -0.4,
                     -0.4, 2.1), 2, 2, byrow = TRUE),
    alpha = c( 0.7,  1.8)
  ),
  list(
    name  = "gen_case_3",
    nu    = 15,
    xi    = c( 0.8,  1.2),
    Omega = matrix(c(2.5, 0.9,
                     0.9, 1.7), 2, 2, byrow = TRUE),
    alpha = c(-1.5,  0.6)
  )
)

## ------------------------------
## 10) run
## ------------------------------
res <- do.call(
  rbind,
  lapply(param_grid, function(pg) {
    ans <- compare_general_skewt(
      xi      = pg$xi,
      Omega   = pg$Omega,
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

res <- res[, c("case", "xi1", "xi2", "omega11", "omega12", "omega22",
               "alpha1", "alpha2", "nu", "k1", "k2", "order",
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
