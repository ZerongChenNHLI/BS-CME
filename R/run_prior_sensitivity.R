## =============================================================================
## Prior sensitivity analysis (Section sim, Table tab:prior).
##
## The main study uses the Dirichlet prior alpha_a = 1/J for every arm (total
## mass 1). This script repeats the selection-stage and end-to-end comparison
## of the co-primary gate (M) and the composite gate with the rule-selected
## weights w*(Theta_conc) under alternative priors, at rho_latent = 0.3, for
## the parallel (S1) and concentrated (S3) configurations:
##
##   weak       Dir(1/J) on every arm, the prior of the paper (mass 1)
##   jeffreys   Dir(1/2) on every arm (mass J/2 = 8)
##   flat       Dir(1) on every arm (mass J = 16; shrinks every marginal to 0.5)
##   ctrl20     doses Dir(1/J); control Dir(20 q_C), q_C the control cell
##              probabilities, i.e. an informative control prior of effective
##              size 20 centred on the truth
##   ctrl50     as ctrl20 with effective size 50
##   ctrl50mis  as ctrl50 but centred on control marginals 0.10 too high
##              (an over-optimistic historical control)
##
## For each prior: n1 at PCS = 0.8 from the PCS curve over n1_grid (common
## stage-1 data for both gates), and at n1 = 100, n_tot = 350 the probability
## of proceeding, PCS, power and E[N], with eta recalibrated under the global
## null for each gate and prior. The prior enters only the selection and
## futility statistics; the confirmatory test is unchanged.
##
##   BSCME_CORES=4 Rscript R/run_prior_sensitivity.R   ->  results/prior_sens.csv
##                                                         results/prior_sens_n1.csv
## About 2 minutes on 4 cores (--quick: a tenth of the replicates).
## =============================================================================

args <- commandArgs(TRUE)
quick <- "--quick" %in% args
cores <- as.integer(Sys.getenv("BSCME_CORES", "1"))
here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else "R"
})
source(file.path(here, "bscme_methods.R"))
source(file.path(here, "scenarios.R"))
out_dir <- file.path(dirname(here), "results")

mc <- list(pcs = 10000, e2e = 10000, eta = 10000)
if (quick) mc <- lapply(mc, function(x) x %/% 10)

RNGkind("L'Ecuyer-CMRG")
set.seed(20261003)
papply <- function(X, FUN) {
  s <- .Random.seed
  seeds <- vector("list", length(X))
  for (i in seq_along(X)) { s <- parallel::nextRNGStream(s); seeds[[i]] <- s }
  f <- function(i) { assign(".Random.seed", seeds[[i]], envir = globalenv()); FUN(X[[i]]) }
  out <- if (cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(seq_along(X), f, mc.cores = cores, mc.set.seed = FALSE)
  } else lapply(seq_along(X), f)
  assign(".Random.seed", parallel::nextRNGStream(s), envir = globalenv())
  out
}
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
write_out <- function(df, name) {
  write.csv(df, file.path(out_dir, name), row.names = FALSE)
  say("wrote", file.path("results", name))
}

## Rule-selected weights from the main run (floored calibration over Theta_conc)
wj <- read.csv(file.path(out_dir, "weights_joint.csv"))
w_conc <- unlist(wj[wj$set == "concordant_floor", paste0("w", 1:4)])
M <- 3; A <- M + 1; J <- 2^K; rho <- 0.3; n1_ref <- 100

## Priors as (M + 1) x J matrices, control last
q_C <- cell_probs_latent(pC, rho = rho)                  # control cell probabilities
q_mis <- cell_probs_latent(pC + 0.10, rho = rho)
prior_mat <- function(dose, ctrl) rbind(matrix(dose, M, J, byrow = TRUE), ctrl)
priors <- list(
  weak      = prior_mat(rep(1 / J, J), rep(1 / J, J)),
  jeffreys  = prior_mat(rep(1 / 2, J), rep(1 / 2, J)),
  flat      = prior_mat(rep(1, J), rep(1, J)),
  ctrl20    = prior_mat(rep(1 / J, J), 20 * q_C),
  ctrl50    = prior_mat(rep(1 / J, J), 50 * q_C),
  ctrl50mis = prior_mat(rep(1 / J, J), 50 * q_mis)
)
prior_mass <- sapply(priors, function(a) paste(round(rowSums(a), 1)[c(1, A)], collapse = "/"))
cfgs <- c("S1", "S3")

## 1. Selection stage: PCS curves and n1 at PCS = 0.8
say("PCS curves by prior")
grid <- expand.grid(prior = names(priors), cfg = cfgs, n1 = n1_grid, stringsAsFactors = FALSE)
pcs <- do.call(rbind, papply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  r <- pcs_selection(pC, configs[[g$cfg]], g$n1, rho, weights = list(C = w_conc),
                     nsim = mc$pcs, alpha0 = priors[[g$prior]])
  data.frame(g, pcs_M = r$pcs[["M"]], pcs_C = r$pcs[["C"]])
}))
n1_tab <- do.call(rbind, lapply(split(pcs, list(pcs$prior, pcs$cfg)), function(d) {
  d <- d[order(d$n1), ]
  data.frame(prior = d$prior[1], cfg = d$cfg[1],
             n1_M = n1_at_target(d$n1, d$pcs_M, pcs_target),
             n1_C = n1_at_target(d$n1, d$pcs_C, pcs_target),
             se_n1_M = n1_se(d$n1, d$pcs_M, pcs_target, mc$pcs),
             se_n1_C = n1_se(d$n1, d$pcs_C, pcs_target, mc$pcs))
}))
n1_tab$ratio <- n1_tab$n1_C / n1_tab$n1_M
write_out(pcs, "prior_sens_pcs_curves.csv")

## 2. End-to-end at n1 = 100, n_tot = 350, eta recalibrated per gate and prior
say("end-to-end by prior")
cells <- expand.grid(prior = names(priors), cfg = cfgs, design = c("M", "C"),
                     stringsAsFactors = FALSE)
e2e <- do.call(rbind, papply(seq_len(nrow(cells)), function(i) {
  g <- cells[i, ]
  gate <- if (g$design == "M") "conjunctive" else "composite"
  w <- if (g$design == "M") NULL else w_conc
  a0 <- priors[[g$prior]]
  eta <- calibrate_eta(pC, M, n1_ref, rho, gate = gate, w = w, target = proceed_null,
                       nsim = mc$eta, alpha0 = a0)
  r <- run_trial(pC, configs[[g$cfg]], n1_ref, n_tot = n_tot, rho = rho, gate = gate,
                 w = w, eta = eta, nsim = mc$e2e, alpha = alpha, alpha0 = a0)
  data.frame(g, n1 = n1_ref, ntot = n_tot, rho = rho, nsim = mc$e2e, eta = eta,
             proceed = r$proceed, pcs = r$pcs, power = r$power, EN = r$EN,
             false_claim = r$false_claim,
             se_pcs = sqrt(r$pcs * (1 - r$pcs) / mc$e2e),
             se_power = sqrt(r$power * (1 - r$power) / mc$e2e))
}))
out <- merge(n1_tab, reshape(e2e[, c("prior", "cfg", "design", "proceed", "pcs", "power", "EN", "se_pcs", "se_power")],
                             idvar = c("prior", "cfg"), timevar = "design", direction = "wide"),
             by = c("prior", "cfg"))
out$mass <- prior_mass[out$prior]
out <- out[order(match(out$cfg, cfgs), match(out$prior, names(priors))), ]
write_out(out, "prior_sens.csv")
say("done", if (quick) "(quick mode)" else "")
