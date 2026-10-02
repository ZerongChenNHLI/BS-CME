## =============================================================================
## BS-CME: Bayesian seamless phase II/III dose-optimisation design with a
## weighted composite selection endpoint and co-primary confirmatory endpoints.
##
## Method code. Base R only (stats, utils). Section and equation numbers refer
## to manuscript/body.tex.
##
##   Data model ....... cell_patterns(), cell_probs_latent(), simulate_counts()
##   Posterior ........ posterior_moments()                      eq. (beta), (moments)
##   Stage-1 gate ..... gate_stats(), select_dose()              eq. (composite), (gate)
##   Final analysis ... wald_p(), simes_p(), inv_normal(),
##                      final_analysis()                         eq. (comb), Sec. final
##   Calibration ...... calibrate_weights(), calibrate_eta(),
##                      ppos()                                   eq. (weights), (ppos)
##   Theory ........... min_normal_var(), n1_ratio_bracket()     eq. (bracket)
##   Trial ............ run_trial()                              Sec. e2e
## =============================================================================

## ---------------------------------------------------------------------------
## Cells and data generation
## ---------------------------------------------------------------------------

#' The J = 2^K response patterns as a J x K 0/1 matrix.
#' Cell j corresponds to y with j = 1 + sum_k y_k 2^(k-1).
cell_patterns <- function(K) {
  as.matrix(expand.grid(rep(list(0:1), K)))[, seq_len(K), drop = FALSE]
}

#' Cell probabilities of a latent Gaussian threshold model.
#'
#' Y_k = 1{Z_k <= qnorm(p_k)}, with Z equicorrelated N(0, 1) with correlation
#' rho. With Z_k = sqrt(rho) U + sqrt(1 - rho) e_k the cells are conditionally
#' independent given U, so each cell probability is a one-dimensional integral.
#'
#' @param p   marginal success probabilities (length K)
#' @param rho latent equicorrelation in [0, 1)
#' @return    length-J vector of cell probabilities, ordered as cell_patterns(K)
cell_probs_latent <- function(p, rho = 0) {
  K <- length(p)
  Y <- cell_patterns(K)
  cut <- qnorm(p)
  if (rho == 0) {
    return(apply(Y, 1, function(y) prod(ifelse(y == 1, p, 1 - p))))
  }
  s <- sqrt(1 - rho)
  pr <- apply(Y, 1, function(y) {
    integrate(function(u) {
      v <- dnorm(u)
      for (k in seq_len(K)) {
        q <- pnorm((cut[k] - sqrt(rho) * u) / s)
        v <- v * if (y[k] == 1) q else 1 - q
      }
      v
    }, -Inf, Inf, rel.tol = 1e-10)$value
  })
  pr / sum(pr)
}

#' Product-moment (phi) correlation between two components implied by the
#' latent model with common marginal p and latent correlation rho.
phi_from_latent <- function(p, rho) {
  pi2 <- cell_probs_latent(c(p, p), rho)
  p11 <- pi2[4]
  (p11 - p^2) / (p * (1 - p))
}

#' Cell-probability matrix for all arms of a configuration.
#'
#' @param pC    control marginals (length K)
#' @param theta M x K matrix of dose effects theta_{a,k} = p_{a,k} - p_{C,k}
#' @param rho   latent equicorrelation
#' @return      (M + 1) x J matrix; rows 1..M are doses, row M + 1 is control
config_cell_probs <- function(pC, theta, rho = 0) {
  P <- rbind(sweep(theta, 2, pC, "+"), pC)
  if (any(P <= 0 | P >= 1)) stop("marginal probabilities must lie in (0, 1)")
  t(apply(P, 1, cell_probs_latent, rho = rho))
}

#' Multinomial cell counts for nsim trials.
#'
#' @param cellp (arms x J) cell probabilities
#' @param n     subjects per arm
#' @param nsim  number of simulated trials
#' @return      array [nsim, arms, J]
simulate_counts <- function(cellp, n, nsim) {
  A <- nrow(cellp); J <- ncol(cellp)
  out <- array(0L, c(nsim, A, J))
  for (a in seq_len(A)) out[, a, ] <- t(rmultinom(nsim, n, cellp[a, ]))
  out
}

## ---------------------------------------------------------------------------
## Dirichlet posterior: exact marginal moments (eq. beta and eq. moments)
## ---------------------------------------------------------------------------

#' Posterior mean and covariance of the K marginals p_{a,k}, for every
#' simulated trial and arm, under pi_a ~ Dir(alpha0) and multinomial counts.
#'
#' Uses aggregation: alpha'(R_k) and alpha'(R_k cap R_k') are sums of
#' posterior Dirichlet parameters over the relevant cells.
#'
#' @param counts array [nsim, arms, J] (or a matrix [arms, J] for one trial)
#' @param alpha0 prior Dirichlet parameter, scalar or length J (default 1/J)
#' @return list(mean = [nsim, arms, K], cov = [nsim, arms, K, K])
posterior_moments <- function(counts, alpha0 = NULL) {
  if (length(dim(counts)) == 2) counts <- array(counts, c(1, dim(counts)))
  d <- dim(counts); nsim <- d[1]; A <- d[2]; J <- d[3]
  K <- as.integer(round(log2(J)))
  Y <- cell_patterns(K)
  if (is.null(alpha0)) alpha0 <- rep(1 / J, J)
  if (length(alpha0) == 1) alpha0 <- rep(alpha0, J)
  post <- sweep(counts, 3, alpha0, "+")                  # alpha'
  flat <- matrix(post, nsim * A, J)                       # rows = (sim, arm)
  a0 <- rowSums(flat)
  aK <- flat %*% Y                                        # alpha'(R_k)
  m <- aK / a0
  cv <- array(0, c(nsim * A, K, K))
  for (k in seq_len(K)) for (l in k:K) {
    akl <- as.vector(flat %*% (Y[, k] * Y[, l]))         # alpha'(R_k cap R_l)
    v <- (a0 * akl - aK[, k] * aK[, l]) / (a0^2 * (a0 + 1))
    cv[, k, l] <- v; cv[, l, k] <- v
  }
  list(mean = array(m, c(nsim, A, K)), cov = array(cv, c(nsim, A, K, K)))
}

## ---------------------------------------------------------------------------
## Stage 1: composite and conjunctive gates (eq. composite, eq. gate)
## ---------------------------------------------------------------------------

#' Interim statistics for every trial and dose.
#'
#' Normal approximation with exact Dirichlet moments, as in Section saving:
#'  * composite:   z_a = (E[Delta_a] - gamma) / sd(Delta_a), Delta_a = S_a - S_C,
#'                 S_a = w' p_a, G_a = Phi(z_a) = Pr(Delta_a > gamma | x1)
#'  * conjunctive: z_{a,k} per endpoint, ranked by min_k z_{a,k};
#'                 gate min_k Pr(p_{a,k} - p_{C,k} > delta_k | x1) = Phi(min_k z_{a,k})
#'
#' @param mom   output of posterior_moments(); last arm is the control
#' @param w     composite weights (length K, sum 1); NULL for conjunctive only
#' @param delta per-endpoint margins (length K), default 0
#' @param gamma composite margin, default sum(w * delta)
#' @return list(zC = [nsim, M] composite z, zM = [nsim, M] min_k z)
gate_stats <- function(mom, w = NULL, delta = NULL, gamma = NULL) {
  d <- dim(mom$mean); nsim <- d[1]; A <- d[2]; K <- d[3]; M <- A - 1
  if (is.null(delta)) delta <- rep(0, K)
  mC <- mom$mean[, A, , drop = FALSE]; vC <- mom$cov[, A, , , drop = FALSE]
  # per-endpoint z
  zk <- array(0, c(nsim, M, K))
  for (a in seq_len(M)) for (k in seq_len(K)) {
    zk[, a, k] <- (mom$mean[, a, k] - mC[, 1, k] - delta[k]) /
      sqrt(mom$cov[, a, k, k] + vC[, 1, k, k])
  }
  zM <- apply(zk, c(1, 2), min)
  zC <- NULL
  if (!is.null(w)) {
    if (is.null(gamma)) gamma <- sum(w * delta)
    sC_mean <- as.vector(matrix(mC, nsim, K) %*% w)
    sC_var <- quad_form(matrix(vC, nsim, K * K), w)
    zC <- matrix(0, nsim, M)
    for (a in seq_len(M)) {
      ma <- as.vector(matrix(mom$mean[, a, ], nsim, K) %*% w)
      va <- quad_form(matrix(mom$cov[, a, , ], nsim, K * K), w)
      zC[, a] <- (ma - sC_mean - gamma) / sqrt(va + sC_var)
    }
  }
  list(zC = zC, zM = zM, zk = zk)
}

## w' Sigma_i w for each row i of a matrix of vectorised K x K covariances.
quad_form <- function(covflat, w) as.vector(covflat %*% as.vector(outer(w, w)))

#' Dose selection and go/no-go (eq. gate).
#'
#' @param z   [nsim, M] ranking statistics (zC or zM from gate_stats)
#' @param eta go/no-go boundary on the probability scale; NULL = always proceed
#' @return    data.frame(sel = selected dose, G = Phi(z) at sel, go = logical)
select_dose <- function(z, eta = NULL) {
  sel <- max.col(z, ties.method = "random")
  G <- pnorm(z[cbind(seq_len(nrow(z)), sel)])
  go <- if (is.null(eta)) rep(TRUE, length(G)) else G >= eta
  data.frame(sel = sel, G = G, go = go)
}

#' Co-primary-optimal dose a° = argmax_a min_k (theta_{a,k} - delta_k).
optimal_dose <- function(theta, delta = rep(0, ncol(theta))) {
  which.max(apply(sweep(theta, 2, delta), 1, min))
}

#' Composite-optimal dose a^C = argmax_a w' theta_a.
composite_optimal_dose <- function(theta, w) which.max(as.vector(theta %*% w))

#' Weights are concordant for a configuration if both criteria pick the same dose.
is_concordant <- function(theta, w, delta = rep(0, ncol(theta))) {
  optimal_dose(theta, delta) == composite_optimal_dose(theta, w)
}

## ---------------------------------------------------------------------------
## Final analysis: closed testing over doses, inverse normal combination,
## intersection-union over endpoints (Section final, eq. comb)
## ---------------------------------------------------------------------------

#' One-sided Wald p-value for H: p_a - p_C <= delta, per endpoint.
#'
#' @param xa,xc   success counts per endpoint (vectors or matrices with K columns)
#' @param na,nc   sample sizes
wald_p <- function(xa, xc, na, nc, delta = 0) {
  pa <- xa / na; pc <- xc / nc
  se <- sqrt(pa * (1 - pa) / na + pc * (1 - pc) / nc)
  z <- (pa - pc - delta) / se
  z[se == 0] <- ifelse((pa - pc - delta)[se == 0] > 0, Inf, -Inf)
  pnorm(z, lower.tail = FALSE)
}

#' Simes p-value for an intersection hypothesis.
simes_p <- function(p) {
  p <- sort(p); m <- length(p)
  min(1, min(m * p / seq_len(m)))
}

#' Bonferroni p-value for an intersection hypothesis.
bonferroni_p <- function(p) min(1, length(p) * min(p))

#' Inverse normal combination function (eq. comb).
#'
#' p-values are clipped away from 0 and 1 so that a stage-1 p of exactly 0
#' (zero standard error with a positive difference) combined with a stage-2 p
#' of exactly 1 does not produce Inf - Inf = NaN.
inv_normal <- function(p1, p2, v1) {
  eps <- 1e-12
  p1 <- pmin(pmax(p1, eps), 1 - eps); p2 <- pmin(pmax(p2, eps), 1 - eps)
  v2 <- sqrt(1 - v1^2)
  pnorm(v1 * qnorm(p1, lower.tail = FALSE) + v2 * qnorm(p2, lower.tail = FALSE),
        lower.tail = FALSE)
}

#' Row-wise Simes / Bonferroni p-values for an nsim x m matrix of p-values.
#'
#' Simes: min_i (m / i) u_(i). Ranks are obtained from column-pair comparisons
#' (m is at most M * K), which is much faster in R than sorting every row; ties
#' receive the larger rank, which reproduces the sorted-order formula.
row_simes <- function(U) {
  m <- ncol(U)
  if (m == 1) return(pmin(1, U[, 1]))
  r <- matrix(0, nrow(U), m)
  for (c1 in seq_len(m)) for (c2 in seq_len(m)) r[, c1] <- r[, c1] + (U[, c2] <= U[, c1])
  pmin(1, do.call(pmin, lapply(seq_len(m), function(c1) m * U[, c1] / r[, c1])))
}
row_bonferroni <- function(U) pmin(1, ncol(U) * do.call(pmin, lapply(seq_len(ncol(U)), function(c1) U[, c1])))

#' Row-wise sort in decreasing order (nsim x m matrix).
row_sort_dec <- function(U) {
  if (ncol(U) == 1) return(U)
  matrix(U[order(row(U), -U)], nrow = nrow(U), byrow = TRUE)
}

#' Co-primary final analysis for nsim trials at once (Section final, Theorem 1).
#'
#' closure = "full" (the design of the manuscript): closed test over the M * K
#' elementary hypotheses H_{a,k}. For an intersection J of dose-endpoint pairs,
#' the stage-1 p-value is the Simes (or Bonferroni) combination of the stage-1
#' p-values of every pair in J, the stage-2 p-value is the combination of the
#' selected dose's stage-2 p-values on the endpoints that J contains for it
#' (1 if there are none), and the two are combined by the inverse-normal rule.
#' H_{a*,k} is rejected iff every J containing (a*, k) is rejected; the claim
#' is made iff H_{a*,k} is rejected for every k. This controls the false-claim
#' rate for every configuration (Theorem 1, Appendix B).
#'
#' The maximum over the 2^{(M-1)K} subsets of the other doses' pairs is taken
#' in closed form: the combination functions are non-decreasing in each
#' argument, so for a given number j of other-dose pairs the largest
#' intersection p-value is attained by the j largest of those p-values. The
#' adjusted p-value for H_{a*,k} is therefore
#'   max over J* (subsets of endpoints containing k) and j in 0..(M-1)K of
#'   C( psi(p1[a*, J*], top_j(others)), psi(p2[J*]) ).
#'
#' closure = "union": the dose-level alternative discussed in the manuscript.
#' The hypotheses are H_a^U = union_k H_{a,k}, each intersection over doses is
#' tested by an intersection-union test whose endpoint-k component combines
#' p1[a*, k] with the dose-level p-values q_a = max_k p1[a, k] of the other
#' doses in the intersection, and the stage-2 p-value is p2[k]. Also valid for
#' every configuration, and more powerful, but it controls a weaker error rate
#' (the false-claim rate only, not the family-wise rate over all H_{a,k}).
#'
#' closure = "endpoint": closes each endpoint separately over doses, using
#' p1[a, k] for the other doses. Controls the false-claim rate only when all
#' not-fully-effective doses share a null endpoint; kept to quantify the
#' inflation in configurations N3 and N4.
#'
#' @param P1    nsim x M x K array of stage-1 p-values (dose a vs stage-1 control)
#' @param P2    nsim x K matrix of stage-2 p-values (selected dose vs stage-2 control)
#' @param sel   length-nsim selected dose
#' @param v1    stage-1 combination weight sqrt(n1 / (n1 + n2_plan))
#' @param alpha one-sided level
#' @param intersection "simes" or "bonferroni"
#' @param closure "full" (default), "union" or "endpoint", see above
#' @return list(padj = nsim x K adjusted p-values for H_{a*,k}, claim = logical nsim)
closed_test <- function(P1, P2, sel, v1, alpha = 0.025,
                        intersection = c("simes", "bonferroni"),
                        closure = c("full", "union", "endpoint")) {
  psi <- switch(match.arg(intersection), simes = row_simes, bonferroni = row_bonferroni)
  closure <- match.arg(closure)
  nsim <- dim(P1)[1]; M <- dim(P1)[2]; K <- dim(P1)[3]
  P2 <- matrix(P2, nsim, K)
  idx <- cbind(seq_len(nsim), sel)
  p1sel <- vapply(seq_len(K), function(k) P1[cbind(idx, k)], numeric(nsim))   # nsim x K
  p1sel <- matrix(p1sel, nsim, K)
  others <- lapply(seq_len(K), function(k) {                                 # per endpoint: nsim x (M-1)
    o <- vapply(seq_len(M), function(a) P1[, a, k], numeric(nsim))
    o <- matrix(o, nsim, M)
    matrix(t(o)[t(col(o) != sel)], nsim, M - 1, byrow = TRUE)
  })
  padj <- matrix(0, nsim, K)
  if (closure == "full") {
    O <- row_sort_dec(do.call(cbind, others))                                # nsim x (M-1)K, decreasing
    Jstar <- unlist(lapply(seq_len(K), function(r) combn(K, r, simplify = FALSE)), recursive = FALSE)
    for (J in Jstar) {
      p2J <- psi(P2[, J, drop = FALSE])
      best <- rep(0, nsim)
      for (j in 0:ncol(O)) {
        U <- if (j == 0) p1sel[, J, drop = FALSE] else cbind(p1sel[, J, drop = FALSE], O[, seq_len(j), drop = FALSE])
        best <- pmax(best, inv_normal(psi(U), p2J, v1))
      }
      for (k in J) padj[, k] <- pmax(padj[, k], best)
    }
  } else {
    q <- row_sort_dec(do.call(pmax, lapply(seq_len(K), function(k) others[[k]])))  # dose-level q_a, decreasing
    for (k in seq_len(K)) {
      Ok <- if (closure == "union") q else row_sort_dec(others[[k]])
      best <- rep(0, nsim)
      for (j in 0:ncol(Ok)) {
        U <- if (j == 0) p1sel[, k, drop = FALSE] else cbind(p1sel[, k], Ok[, seq_len(j), drop = FALSE])
        best <- pmax(best, inv_normal(psi(U), P2[, k], v1))
      }
      padj[, k] <- best
    }
  }
  list(padj = padj, claim = rowSums(padj <= alpha) == K)
}

#' Co-primary final analysis for one trial: wrapper around closed_test().
#'
#' @param p1    M x K matrix of stage-1 p-values (dose a vs stage-1 control)
#' @param p2    length-K stage-2 p-values (selected dose vs stage-2 control)
#' @param sel   selected dose
#' @return list(reject = length-K logical, claim = all(reject), padj = closed-test p)
final_analysis <- function(p1, p2, sel, v1, alpha = 0.025,
                           intersection = c("simes", "bonferroni"),
                           closure = c("full", "union", "endpoint")) {
  M <- nrow(p1); K <- ncol(p1)
  ct <- closed_test(array(p1, c(1, M, K)), matrix(p2, 1, K), sel, v1, alpha,
                    intersection, closure)
  padj <- as.vector(ct$padj)
  list(reject = padj <= alpha, claim = all(padj <= alpha), padj = padj)
}

#' Reference implementation of the full closure by explicit enumeration of all
#' 2^(MK) - 1 intersections, for checking closed_test() on small problems.
final_analysis_enum <- function(p1, p2, sel, v1, alpha = 0.025,
                                intersection = c("simes", "bonferroni")) {
  f <- match.fun(paste0(match.arg(intersection), "_p"))
  M <- nrow(p1); K <- ncol(p1); L <- M * K
  pa <- rep(seq_len(M), times = K); pk <- rep(seq_len(K), each = M)      # column a + M(k-1)
  p1v <- as.vector(p1)
  B <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), L)))[-1, , drop = FALSE]
  cp <- apply(B, 1, function(b) {
    s1 <- f(p1v[b]); ks <- pk[b & pa == sel]
    s2 <- if (length(ks)) f(p2[ks]) else 1
    inv_normal(s1, s2, v1)
  })
  padj <- vapply(seq_len(K), function(k) max(cp[B[, sel + M * (k - 1)]]), numeric(1))
  list(reject = padj <= alpha, claim = all(padj <= alpha), padj = padj)
}

## ---------------------------------------------------------------------------
## Design calibration (Section calib)
## ---------------------------------------------------------------------------

#' All weight vectors on the K-simplex with the given grid step.
simplex_grid <- function(K, step = 0.05) {
  n <- round(1 / step)
  g <- as.matrix(expand.grid(rep(list(0:n), K - 1)))
  g <- g[rowSums(g) <= n, , drop = FALSE]
  W <- cbind(g, n - rowSums(g)) / n
  unname(W)
}

#' Calibrated weights (eq. weights): maximise mean PCS over a scenario set.
#'
#' Common random numbers: the stage-1 data for each scenario are simulated
#' once and every weight vector is evaluated on them.
#'
#' @param scenarios list of M x K theta matrices
#' @param pC        control marginals
#' @param n1        stage-1 size per arm
#' @param rho       latent correlation
#' @param W         candidate weights (rows); default simplex grid, step 0.05
#' @param constraint optional function(w) -> TRUE/FALSE encoding the set cal{W}
#' @return list(w = best weights, pcs = its mean PCS, all = PCS per candidate)
calibrate_weights <- function(scenarios, pC, n1, rho = 0, W = NULL,
                              nsim = 2000, delta = NULL, alpha0 = NULL,
                              constraint = NULL) {
  K <- length(pC)
  if (is.null(W)) W <- simplex_grid(K)
  if (!is.null(constraint)) W <- W[apply(W, 1, constraint), , drop = FALSE]
  if (is.null(delta)) delta <- rep(0, K)
  pcs <- matrix(0, nrow(W), length(scenarios))
  for (s in seq_along(scenarios)) {
    theta <- scenarios[[s]]
    M <- nrow(theta)
    best <- optimal_dose(theta, delta)
    cnt <- simulate_counts(config_cell_probs(pC, theta, rho), n1, nsim)
    mom <- posterior_moments(cnt, alpha0)
    mC <- matrix(mom$mean[, M + 1, ], nsim, K)
    vC <- matrix(mom$cov[, M + 1, , ], nsim, K * K)
    WW <- t(apply(W, 1, function(w) as.vector(outer(w, w))))  # nW x K^2
    gam <- as.vector(W %*% delta)
    z <- array(0, c(nsim, M, nrow(W)))
    for (a in seq_len(M)) {
      dm <- matrix(mom$mean[, a, ], nsim, K) - mC               # nsim x K
      vv <- matrix(mom$cov[, a, , ], nsim, K * K) + vC           # nsim x K^2
      z[, a, ] <- (dm %*% t(W) - rep(gam, each = nsim)) / sqrt(vv %*% t(WW))
    }
    for (i in seq_len(nrow(W))) {
      pcs[i, s] <- mean(max.col(z[, , i], ties.method = "random") == best)
    }
  }
  m <- rowMeans(pcs)
  list(w = W[which.max(m), ], pcs = max(m), all = data.frame(W, pcs = m))
}

#' Futility boundary calibrated so that Pr(proceed | global null) = target.
#'
#' @param gate "composite" or "conjunctive"
calibrate_eta <- function(pC, M, n1, rho = 0, gate = c("composite", "conjunctive"),
                          w = NULL, target = 0.20, nsim = 10000,
                          delta = NULL, alpha0 = NULL) {
  gate <- match.arg(gate)
  K <- length(pC)
  theta0 <- matrix(0, M, K)
  cnt <- simulate_counts(config_cell_probs(pC, theta0, rho), n1, nsim)
  gs <- gate_stats(posterior_moments(cnt, alpha0), w = w, delta = delta)
  z <- if (gate == "composite") gs$zC else gs$zM
  G <- pnorm(apply(z, 1, max))
  unname(quantile(G, 1 - target, type = 7))
}

#' Predictive probability of co-primary success for stage-2 size n2 (eq. ppos).
#'
#' Monte Carlo over Dirichlet posteriors of the selected dose and control,
#' simulating stage-2 data and running the final analysis.
#'
#' @param x1   (M + 1) x J stage-1 counts of a single trial
#' @param sel  selected dose
ppos <- function(x1, sel, n1, n2, n2_plan = n2, ndraw = 2000, alpha = 0.025,
                 delta = NULL, alpha0 = NULL) {
  J <- ncol(x1); K <- as.integer(round(log2(J))); A <- nrow(x1)
  Y <- cell_patterns(K)
  if (is.null(alpha0)) alpha0 <- rep(1 / J, J)
  if (is.null(delta)) delta <- rep(0, K)
  succ1 <- x1 %*% Y
  p1 <- wald_p(succ1[-A, , drop = FALSE],
               matrix(succ1[A, ], A - 1, K, byrow = TRUE), n1, n1,
               matrix(delta, A - 1, K, byrow = TRUE))
  v1 <- sqrt(n1 / (n1 + n2_plan))
  rdir <- function(a) { g <- rgamma(length(a), a); g / sum(g) }
  hits <- logical(ndraw)
  for (i in seq_len(ndraw)) {
    xa <- rmultinom(1, n2, rdir(alpha0 + x1[sel, ]))
    xc <- rmultinom(1, n2, rdir(alpha0 + x1[A, ]))
    p2 <- wald_p(drop(t(xa) %*% Y), drop(t(xc) %*% Y), n2, n2, delta)
    hits[i] <- final_analysis(p1, p2, sel, v1, alpha)$claim
  }
  mean(hits)
}

#' Smallest n2 in [n_min, n_max] with PPoS >= 1 - beta (eq. ppos).
reestimate_n2 <- function(x1, sel, n1, n_min, n_max, n2_plan, beta = 0.2,
                          step = 10, ...) {
  for (n2 in seq(n_min, n_max, by = step)) {
    if (ppos(x1, sel, n1, n2, n2_plan = n2_plan, ...) >= 1 - beta) return(n2)
  }
  n_max
}

## ---------------------------------------------------------------------------
## Second-moment theory (eq. bracket)
## ---------------------------------------------------------------------------

#' Variance of the minimum of K independent standard normals.
min_normal_var <- function(K) {
  f <- function(x, r) x^r * K * dnorm(x) * pnorm(x, lower.tail = FALSE)^(K - 1)
  m1 <- integrate(f, -Inf, Inf, r = 1)$value
  m2 <- integrate(f, -Inf, Inf, r = 2)$value
  m2 - m1^2
}

#' Bracket for n1^C / n1^M under parallel effects (eq. bracket).
#' @param phi product-moment correlation of the binary components
n1_ratio_bracket <- function(K, phi) {
  num <- phi + (1 - phi) / K
  c(lower = num, upper = num / (phi + (1 - phi) * min_normal_var(K)))
}

## ---------------------------------------------------------------------------
## Selection-stage PCS and end-to-end trial simulation
## ---------------------------------------------------------------------------

#' Probability of correct selection at the interim, selection only.
#'
#' All gates are evaluated on the same simulated stage-1 data (common random
#' numbers), so differences between gates within a call are paired.
#'
#' @param weights named list of weight vectors for composite gates
#' @return list(pcs = named vector of PCS for "M" (conjunctive) and each
#'   composite weight, regret = named vector of the expected regret
#'   E[min_k(theta_{a°,k} - delta_k) - min_k(theta_{a*,k} - delta_k)])
pcs_selection <- function(pC, theta, n1, rho = 0, weights = list(),
                          nsim = 4000, delta = NULL, alpha0 = NULL) {
  K <- length(pC)
  if (is.null(delta)) delta <- rep(0, K)
  crit <- apply(sweep(theta, 2, delta), 1, min)          # co-primary criterion per dose
  best <- which.max(crit)
  cnt <- simulate_counts(config_cell_probs(pC, theta, rho), n1, nsim)
  mom <- posterior_moments(cnt, alpha0)
  sel <- list(M = select_dose(gate_stats(mom, delta = delta)$zM)$sel)
  for (nm in names(weights)) {
    sel[[nm]] <- select_dose(gate_stats(mom, w = weights[[nm]], delta = delta)$zC)$sel
  }
  list(pcs = vapply(sel, function(s) mean(s == best), numeric(1)),
       regret = vapply(sel, function(s) mean(crit[best] - crit[s]), numeric(1)))
}

#' n1 at which a PCS curve first reaches the target, by linear interpolation.
#' Returns Inf if the target is not reached on the grid.
n1_at_target <- function(n1, pcs, target = 0.8) {
  i <- which(pcs >= target)
  if (!length(i)) return(Inf)
  i <- i[1]
  if (i == 1) return(n1[1])
  n1[i - 1] + (target - pcs[i - 1]) * (n1[i] - n1[i - 1]) / (pcs[i] - pcs[i - 1])
}

#' End-to-end simulation of the seamless design (Section e2e).
#'
#' Stage 1 randomises n1 to each of M doses and control; the gate selects a
#' dose and decides go/no-go; stage 2 randomises n2 = n_tot - n1 to the
#' selected dose and control; the final analysis uses Simes closure over
#' doses, inverse normal combination with v1 = sqrt(n1 / n_tot) and the
#' intersection-union claim over endpoints.
#'
#' @param gate "composite" or "conjunctive"
#' @param eta  futility boundary (from calibrate_eta)
#' @param others passed to final_analysis(): "union" (the manuscript's closed
#'   test) or "endpoint" (per-endpoint closure, for the N3 comparison only)
#' @return list of operating characteristics:
#'   proceed      probability of proceeding to stage 2
#'   pcs          probability that the selected dose is the co-primary-optimal dose
#'   power        probability of a correct claim: the claim is made and the
#'                selected dose is effective on every endpoint
#'   claim_opt    probability of a claim on the co-primary-optimal dose
#'   false_claim  probability of a claim on a dose with some H_{a,k} true
#'   EN           expected total sample size
#'   trials       per-trial data
run_trial <- function(pC, theta, n1, n_tot = 350, rho = 0,
                      gate = c("composite", "conjunctive"), w = NULL, eta = 0,
                      nsim = 10000, alpha = 0.025, delta = NULL, alpha0 = NULL,
                      intersection = "simes", closure = "full") {
  gate <- match.arg(gate)
  K <- length(pC); M <- nrow(theta); A <- M + 1
  if (is.null(delta)) delta <- rep(0, K)
  n2 <- n_tot - n1
  v1 <- sqrt(n1 / n_tot)
  cellp <- config_cell_probs(pC, theta, rho)
  Y <- cell_patterns(K)

  # Stage 1
  x1 <- simulate_counts(cellp, n1, nsim)
  gs <- gate_stats(posterior_moments(x1, alpha0), w = w, delta = delta)
  sd <- select_dose(if (gate == "composite") gs$zC else gs$zM, eta)

  # Stage-1 per-endpoint p-values, every dose vs stage-1 control
  s1 <- array(0, c(nsim, A, K))
  for (a in seq_len(A)) s1[, a, ] <- matrix(x1[, a, ], nsim) %*% Y
  p1 <- array(0, c(nsim, M, K))
  for (a in seq_len(M)) {
    p1[, a, ] <- wald_p(s1[, a, ], s1[, A, ], n1, n1,
                        matrix(delta, nsim, K, byrow = TRUE))
  }

  # Stage 2, for the trials that proceed: counts simulated by selected dose,
  # per-endpoint p-values against the stage-2 control, closed test on all
  # proceeding trials at once.
  claim <- logical(nsim)
  null_true <- sweep(theta, 2, delta) <= 0              # H_{a,k} true
  false_claim <- logical(nsim)
  go <- which(sd$go)
  if (length(go)) {
    p2 <- matrix(0, length(go), K)
    for (a in unique(sd$sel[go])) {
      i <- which(sd$sel[go] == a)
      xa <- t(rmultinom(length(i), n2, cellp[a, ])); xc <- t(rmultinom(length(i), n2, cellp[A, ]))
      p2[i, ] <- wald_p(xa %*% Y, xc %*% Y, n2, n2, matrix(delta, length(i), K, byrow = TRUE))
    }
    ct <- closed_test(p1[go, , , drop = FALSE], p2, sd$sel[go], v1, alpha, intersection, closure)
    claim[go] <- ct$claim
    false_claim[go] <- ct$claim & apply(null_true, 1, any)[sd$sel[go]]
  }
  best <- optimal_dose(theta, delta)
  list(proceed = mean(sd$go),
       pcs = mean(sd$sel == best),
       power = mean(claim & !false_claim),
       claim_opt = mean(claim & sd$sel == best),
       false_claim = mean(false_claim),
       EN = A * n1 + 2 * n2 * mean(sd$go),
       trials = data.frame(sel = sd$sel, G = sd$G, go = sd$go, claim = claim,
                           false_claim = false_claim))
}
