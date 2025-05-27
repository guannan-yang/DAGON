# -------------------------------------------------------------------------
#  DRAGON   (Dependence Relationship And Grouping Of Networks)   –  R PORT
#  -----------------------------------------------------------------------
#  • Maintains the original NetZooR pipeline for penalty selection and
#    shrunken covariance / precision / partial‑correlation generation.
#  • Adds a **faithful translation** of the Python kappa–based p‑value
#    system, including:  estimate_kappa(), estimate_kappa_dragon(), and
#    estimate_p_values_dragon().
#  • dragon() now supports   pval = TRUE   to return   p_raw / p_adj
#    using the analytic / estimated kappa null distribution.
#  -----------------------------------------------------------------------
#  CRAN dependencies: MASS  (simulation), stats (base), hypergeo (2F1).
#  Install with:   install.packages(c("MASS","hypergeo"))
#  -----------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required – install.packages('MASS')")
  if (!requireNamespace("hypergeo", quietly = TRUE))
    stop("Package 'hypergeo' (Gauss hypergeometric 2F1) is required – install.packages('hypergeo')")
})

# =============== 1.  CORE HELPERS  (unchanged original) ==================
#  Note: we assume the user has already sourced the original NetZooR
#  code that defines: VarS, EsqS, estimatePenaltyParameters, 
#                    get_shrunken_covariance_dragon, get_precision_matrix_dragon,
#                    get_partial_correlation_dragon.
if (!exists("estimatePenaltyParameters"))
  stop("Please source the original NetZooR 'DRAGON.R' before this file.")

# ---------------- Monte‑Carlo null (used internally by kappa est.) -------
MC_estimate <- function(n, p1, p2, lambdas, seed = 1) {
  set.seed(seed)
  X1 <- MASS::mvrnorm(n, mu = rep(0, p1), Sigma = diag(p1))
  X2 <- MASS::mvrnorm(n, mu = rep(0, p2), Sigma = diag(p2))
  get_partial_correlation_dragon(X1, X2, lambdas)
}

# ================= 2.  SINGLE‑LAYER kappa  ===============================
logli_term <- function(rho, lambda0) 0.5 * log1p(- (rho / (1 - lambda0))^2)

estimate_kappa <- function(n, p, lambda0, seed = 1) {
  set.seed(seed)
  X  <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  r  <- get_partial_correlation_from_precision(
    get_precision_matrix_dragon(X, matrix(,0,0), lambdas = lambda0)[[1]])
  r  <- r[upper.tri(r)]
  term <- sum(logli_term(r, lambda0))
  m    <- length(r)
  Dlogli <- function(k) 0.5 * m * (digamma(k / 2) - digamma((k - 1) / 2)) + term
  uniroot(Dlogli, interval = c(1.001, 1000 * n))$root
}

# ================= 3.  TWO‑LAYER kappa  ==================================
logli_term_lam <- function(rho, lam) 0.5 * log1p(- (rho / lam)^2)

estimate_kappa_dragon <- function(n, p1, p2, lambdas, seed = 1,
                                  simultaneous = FALSE) {
  set.seed(seed)
  X1 <- MASS::mvrnorm(n, mu = rep(0, p1), Sigma = diag(p1))
  X2 <- MASS::mvrnorm(n, mu = rep(0, p2), Sigma = diag(p2))
  r_sim <- get_partial_correlation_dragon(X1, X2, lambdas)
  
  # split vectors ---------------------------------------------------------
  r11 <- r_sim[1:p1, 1:p1][upper.tri(diag(p1))]
  r22 <- r_sim[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)][upper.tri(diag(p2))]
  r12 <- as.vector(r_sim[1:p1, (p1 + 1):(p1 + p2)])
  
  lam1 <- lambdas[1]; lam2 <- lambdas[2]
  k11 <- k22 <- k12 <- NA
  
  if (lam1 == 1 || lam2 == 1) simultaneous <- FALSE  # match Python guard
  
  if (!simultaneous) {
    if (lam1 < 1) {
      term <- sum(logli_term_lam(r11, 1 - lam1))
      Dlogli <- function(k) 0.25 * p1 * (p1 - 1) * (digamma(k / 2) -
                                                      digamma((k - 1) / 2)) + term
      k11 <- uniroot(Dlogli, c(1.001, 1000 * n))$root
    }
    if (lam2 < 1) {
      term <- sum(logli_term_lam(r22, 1 - lam2))
      Dlogli <- function(k) 0.25 * p2 * (p2 - 1) * (digamma(k / 2) -
                                                      digamma((k - 1) / 2)) + term
      k22 <- uniroot(Dlogli, c(1.001, 1000 * n))$root
    }
    if (lam1 < 1 && lam2 < 1) {
      term <- sum(logli_term_lam(r12, sqrt((1 - lam1) * (1 - lam2))))
      Dlogli <- function(k) 0.25 * 2 * p1 * p2 * (digamma(k / 2) -
                                                    digamma((k - 1) / 2)) + term
      k12 <- uniroot(Dlogli, c(1.001, 1000 * n))$root
    }
  } else {
    term <- sum(logli_term_lam(r11, 1 - lam1)) +
      sum(logli_term_lam(r22, 1 - lam2)) +
      sum(logli_term_lam(r12, sqrt((1 - lam1) * (1 - lam2))))
    coeff <- 0.25 * (p1 + p2) * (p1 + p2 - 1)
    Dlogli <- function(k) coeff * (digamma(k / 2) - digamma((k - 1) / 2)) + term
    k_all <- uniroot(Dlogli, c(1.001, 1000 * n))$root
    k11 <- k22 <- k12 <- k_all
  }
  
  c(kappa11 = k11, kappa22 = k22, kappa12 = k12)
}

# ================= 4.  DRAGON p‑values  ==================================
#  hypergeo::hypergeo(a,b,c,z) gives 2F1  (complex‑safe; we use Re)

pfunc <- function(x, l, k) as.numeric(x * Re(hypergeo::hypergeo(0.5, (3 - k) / 2,
                                                                1.5, (x / l)^2)))

denominator <- function(l, k) if (!is.na(k)) beta(0.5, (k - 1) / 2) * l else NA

pval_scalar <- function(r, l, k, d) {
  res <- 1 - 2 * pfunc(abs(r), l, k) / d
  if (abs(res) < 1e-8) {
    f <- function(r0) 2 * (1 - (r0^2) / l^2)^(k / 2 - 1.5) / d
    res <- 2 * integrate(f, lower = abs(r), upper = l)$value  # numeric fallback
  }
  res
}

pval_vec <- function(r_vec, l, k, d) vapply(r_vec, pval_scalar, numeric(1), l = l, k = k, d = d)

estimate_p_values_dragon <- function(r, n, p1, p2, lambdas,
                                     kappa = "estimate", seed = 1,
                                     simultaneous = FALSE) {
  p_total <- p1 + p2
  lam1 <- lambdas[1]; lam2 <- lambdas[2]
  
  # split observed r ------------------------------------------------------
  r11 <- r[1:p1, 1:p1][upper.tri(diag(p1))]
  r22 <- r[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)][upper.tri(diag(p2))]
  r12 <- as.vector(r[1:p1, (p1 + 1):(p1 + p2)])
  
  if (identical(kappa, "estimate")) {
    if (lam1 == 1 && lam2 == 1) {
      p_raw <- matrix(1, p_total, p_total)
      return(list(p_raw = p_raw, p_adj = p_raw))
    }
    kap_samples <- replicate(10, estimate_kappa_dragon(n, p1, p2, lambdas,
                                                       seed = seed + sample.int(1e6,1),
                                                       simultaneous = simultaneous))
    kappa_est <- rowMeans(kap_samples)
  } else if (identical(kappa, "analytic")) {
    kappa_scalar <- n - 1 - (p_total - 2)
    kappa_est <- rep(kappa_scalar, 3)
  } else {
    kappa_est <- kappa  # assume numeric vector length 3
  }
  
  # denoms ---------------------------------------------------------------
  d11 <- denominator(1 - lam1, kappa_est[1])
  d22 <- denominator(1 - lam2, kappa_est[2])
  d12 <- denominator(sqrt((1 - lam1) * (1 - lam2)), kappa_est[3])
  
  n1 <- p1 * (p1 - 1) / 2
  n2 <- p2 * (p2 - 1) / 2
  n12 <- p1 * p2
  
  if (lam1 < 1) {
    p11 <- pval_vec(r11, 1 - lam1, kappa_est[1], d11)
  } else p11 <- rep(1, n1)
  
  if (lam2 < 1) {
    p22 <- pval_vec(r22, 1 - lam2, kappa_est[2], d22)
  } else p22 <- rep(1, n2)
  
  if (lam1 < 1 && lam2 < 1) {
    l12 <- sqrt((1 - lam1) * (1 - lam2))
    p12 <- pval_vec(r12, l12, kappa_est[3], d12)
  } else p12 <- rep(1, n12)
  
  # assemble -------------------------------------------------------------
  p_raw <- matrix(0, p_total, p_total)
  p_raw[1:p1, 1:p1][upper.tri(diag(p1))] <- p11
  p_raw[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)][upper.tri(diag(p2))] <- p22
  p_raw[1:p1, (p1 + 1):(p1 + p2)] <- matrix(p12, nrow = p1)
  p_raw <- p_raw + t(p_raw)
  
  # adjustment -----------------------------------------------------------
  adj11 <- p.adjust(p11, method = "BH")
  adj22 <- p.adjust(p22, method = "BH")
  adj12 <- p.adjust(p12, method = "BH")
  p_adj <- matrix(0, p_total, p_total)
  p_adj[1:p1, 1:p1][upper.tri(diag(p1))] <- adj11
  p_adj[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)][upper.tri(diag(p2))] <- adj22
  p_adj[1:p1, (p1 + 1):(p1 + p2)] <- matrix(adj12, nrow = p1)
  p_adj <- p_adj + t(p_adj)
  
  list(p_raw = p_raw, p_adj = p_adj,
       kappa = kappa_est)
}

# ================= 5.  WRAPPER ==========================================
#' DRAGON partial‑correlation network with kappa‑based p‑values
#' @param layer1,layer2 paired omics matrices (n × p1 / p2)
#' @param pval logical – compute p‑values? (default FALSE)
#' @param kappa choice: 'estimate' (default), 'analytic', or numeric vector
#' @param simultaneous logical – fit kappa parameters jointly? (default FALSE)
#' @param seed RNG seed
#' @param adj_method multiple‑testing adjustment (default "BH")
#' @return list   cov, prec, ggm, lambdas, gammas, risk_grid, p_raw, p_adj
#' @export

dragon <- function(layer1, layer2, pval = T,
                   kappa = "estimate", simultaneous = FALSE,
                   seed = 1, adj_method = "BH", fdr = 0.1, verbose = FALSE) {
  if (verbose) message("[DRAGON] Estimating shrinkage parameters …")
  params <- estimatePenaltyParameters(layer1, layer2)
  lambdas <- params$lambdas
  
  if (verbose) message("[DRAGON] Building shrunken matrices …")
  cov_mat  <- get_shrunken_covariance_dragon(layer1, layer2, lambdas)
  prec_mat <- get_precision_matrix_dragon(layer1, layer2, lambdas)
  ggm_mat  <- get_partial_correlation_dragon(layer1, layer2, lambdas)
  
  out <- list(cov        = cov_mat,
              prec       = prec_mat,
              ggm        = ggm_mat,
              lambdas    = lambdas,
              gammas     = params$gammas,
              risk_grid  = params$risk_grid)
  
  if (pval) {
    if (verbose) message("[DRAGON] Calculating p‑values …")
    n  <- nrow(layer1)
    p1 <- ncol(layer1)
    p2 <- ncol(layer2)
    
    pv <- estimate_p_values_dragon(ggm_mat, n, p1, p2, lambdas,
                                   kappa = kappa, seed = seed,
                                   simultaneous = simultaneous)
    out$p_raw <- pv$p_raw
    out$p_adj <- pv$p_adj
    out$kappa <- pv$kappa
    
    selection <- matrix(FALSE, nrow = p1 + p2, ncol = p1 + p2)
    selection[which(out$p_adj < fdr, arr.ind = TRUE)] <- TRUE
    diag(selection) <- FALSE  # no self‑connections
    out$selection <- selection
  }
  
  class(out) <- "dragon"
  out
}

