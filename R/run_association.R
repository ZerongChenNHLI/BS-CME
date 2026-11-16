## =============================================================================
## Agreement between the interim composite and the co-primary criterion
## (Table tab:assoc). The endpoint switch rests on the composite ranking the
## doses as the co-primary criterion would (concordance, Section sec:saving).
## This script measures how closely the two agree on finite stage-1 data.
##
## For each configuration S1-S4, rho_lat in {0, 0.3, 0.6} and weight vector
## (uniform, rule-selected w*(Theta_conc)), n1 = 100 per arm, n_tot = 350,
## futility switched off so that every trial proceeds, 10,000 trials:
##   r_z      Pearson correlation between the composite z-statistic zC and the
##            conjunctive statistic zM = min_k z_k over all (trial, dose) pairs
##   tau      mean Kendall tau between the within-trial dose rankings by zC and zM
##   agree    Pr(the composite and the co-primary gate select the same dose)
##   agree_opt Pr(both select the co-primary-optimal dose)
##   r_G      correlation, over trials, between G = Pr(composite effect > 0 | x1)
##            and GM = min_k Pr(p_{a*,k} > p_{C,k} | x1) at the composite-selected dose
##   auc_G    AUC of G for the eventual co-primary claim on the selected dose
##   auc_GM   AUC of GM for the same claim
##   claim    probability of that claim (any dose; power + false claim)
## and, for the co-primary gate on the same configurations, the AUC of its own
## statistic GM at its own selected dose (design = "M", auc_GM).
##
##   Rscript R/run_association.R        -> results/association.csv (about 2 minutes)
## Own random-number seed; nothing else is changed.
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
nsim <- if (quick) 1000 else 10000

RNGkind("L'Ecuyer-CMRG")
set.seed(20261005)
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

wj <- read.csv(file.path(out_dir, "weights_joint.csv"))
w_conc <- unlist(wj[wj$set == "concordant_floor", paste0("w", 1:4)])
weights <- list(Cu = w_uniform, Cj_conc = w_conc)

auc <- function(score, y) {               # Mann-Whitney AUC, ties counted 1/2
  y <- as.logical(y)
  if (!any(y) || all(y)) return(NA_real_)
  r <- rank(score)
  (sum(r[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
}
kendall_rows <- function(X, Y) {          # mean Kendall tau over rows (M = 3 doses)
  M <- ncol(X); pr <- combn(M, 2)
  conc <- rowSums(sapply(seq_len(ncol(pr)), function(j)
    sign(X[, pr[1, j]] - X[, pr[2, j]]) * sign(Y[, pr[1, j]] - Y[, pr[2, j]])))
  mean(conc / ncol(pr))
}

cells <- expand.grid(cfg = names(configs), rho = rho_latent, design = c(names(weights), "M"),
                     stringsAsFactors = FALSE)
say("simulating", nrow(cells), "cells")
res <- do.call(rbind, papply(seq_len(nrow(cells)), function(i) {
  g <- cells[i, ]
  th <- configs[[g$cfg]]
  best <- optimal_dose(th)
  if (g$design == "M") {
    r <- run_trial(pC, th, n1 = 100, n_tot = n_tot, rho = g$rho, gate = "conjunctive",
                   eta = 0, nsim = nsim, alpha = alpha)
    return(data.frame(cfg = g$cfg, rho = g$rho, design = "M", n1 = 100, nsim = nsim,
                      r_z = NA, tau = NA, agree = NA, agree_opt = NA, r_G = NA,
                      auc_G = NA, auc_GM = auc(r$trials$G, r$trials$claim),
                      pcs = r$pcs, claim = mean(r$trials$claim)))
  }
  r <- run_trial(pC, th, n1 = 100, n_tot = n_tot, rho = g$rho, gate = "composite",
                 w = weights[[g$design]], eta = 0, nsim = nsim, alpha = alpha)
  zC <- r$stats$zC; zM <- r$stats$zM
  selC <- max.col(zC, ties.method = "first"); selM <- max.col(zM, ties.method = "first")
  data.frame(cfg = g$cfg, rho = g$rho, design = g$design, n1 = 100, nsim = nsim,
             r_z = cor(as.vector(zC), as.vector(zM)),
             tau = kendall_rows(zC, zM),
             agree = mean(selC == selM),
             agree_opt = mean(selC == best & selM == best),
             r_G = cor(r$trials$G, r$trials$GM),
             auc_G = auc(r$trials$G, r$trials$claim),
             auc_GM = auc(r$trials$GM, r$trials$claim),
             pcs = r$pcs, claim = mean(r$trials$claim))
}))
res <- res[order(match(res$cfg, names(configs)), res$rho, match(res$design, c("Cu", "Cj_conc", "M"))), ]
write.csv(res, file.path(out_dir, "association.csv"), row.names = FALSE)
say("wrote results/association.csv", if (quick) "(quick mode)" else "")
