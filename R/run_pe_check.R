## =============================================================================
## Composite gate with a per-endpoint futility check (Section e2e; extra rows in
## Tables e2eA, e2eN and e2eD; N0/N1 columns of Table e2eD; Figure fig:frontier).
##
## The composite gate ranks and gates on the weighted composite as before. In
## addition the trial stops at the interim if the selected dose has
##   min_k Pr(p_{a*,k} - p_{C,k} > 0 | stage-1 data) < eta_pe,
## with eta_pe = 0.5: a stop whenever the posterior mean difference is negative
## on some endpoint. The composite boundary eta is then recalibrated so that the
## two rules together proceed with probability 0.20 under the global null.
## Theorem 1 covers any futility rule, so the confirmatory test is unchanged.
##
## Designs added: Cu_pe (uniform weights + check), Cj_conc_pe (rule-selected
## weights w*(Theta_conc) + check). Reference designs M, Cu, Cj_conc, Cc_S3 are
## re-created here only for block N.
##
## Block A  n1 = 100, n_tot = 350, rho in {0, 0.3, 0.6}, S1-S4 and N0-N4
##          -> results/e2e_pe_runA.csv (same columns as e2e_runA.csv)
## Block D  (n1, n_tot) grid at rho = 0.3 for S1 (Cu_pe, Cj_conc_pe) and S3
##          (Cj_conc_pe), optimised for target power 0.80/0.85/0.90 against the
##          co-primary gate of e2e_runD.csv
##          -> results/e2e_pe_runD.csv, results/e2e_pe_runD_matched.csv
## Block N  E[N] under N0 and N1 at every design's run-D optimum (n1, n_tot),
##          for the designs of Table e2eD and the two new ones
##          -> results/e2e_runD_nulls.csv
## Redraws results/figures/sim_e2e_fig.pdf with the two new curves.
##
##   BSCME_CORES=4 Rscript R/run_pe_check.R        (about 5 minutes; --quick: a tenth)
## Requires results/weights_joint.csv, weights_calibrated.csv, e2e_runD.csv and
## e2e_runD_matched.csv
## from run_simulations.R. Own random-number seed; nothing else is changed.
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
fig_dir <- file.path(out_dir, "figures")

mc <- list(e2e = 10000, e2e_null = 40000, eta = 10000, nulls = 20000)
if (quick) mc <- lapply(mc, function(x) x %/% 10)

RNGkind("L'Ecuyer-CMRG")
set.seed(20261004)
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
rd <- function(name) read.csv(file.path(out_dir, name))

eta_pe <- 0.5
wj <- rd("weights_joint.csv")
w_conc <- unlist(wj[wj$set == "concordant_floor", paste0("w", 1:4)])
wcal <- rd("weights_calibrated.csv")
w_star <- function(cfg, rho) unlist(wcal[wcal$cfg == cfg & wcal$rho == rho, paste0("w", 1:4)])
designs <- list(
  M          = list(gate = "conjunctive", w = NULL,              eta_pe = NULL),
  Cu         = list(gate = "composite",   w = w_uniform,         eta_pe = NULL),
  Cj_conc    = list(gate = "composite",   w = w_conc,            eta_pe = NULL),
  Cc_S3      = list(gate = "composite",   w = w_star("S3", 0.3), eta_pe = NULL),
  Cu_pe      = list(gate = "composite",   w = w_uniform,         eta_pe = eta_pe),
  Cj_conc_pe = list(gate = "composite",   w = w_conc,            eta_pe = eta_pe)
)
new_designs <- c("Cu_pe", "Cj_conc_pe")
all_configs <- c(configs, null_configs)

eta_cache <- new.env()
get_eta <- function(design, n1, rho = 0.3) {
  key <- paste(design, n1, rho)
  if (is.null(eta_cache[[key]])) {
    d <- designs[[design]]
    eta_cache[[key]] <- calibrate_eta(pC, M = 3, n1 = n1, rho = rho, gate = d$gate, w = d$w,
                                      target = proceed_null, nsim = mc$eta, eta_pe = d$eta_pe)
  }
  eta_cache[[key]]
}
e2e_cell <- function(cfg, design, n1, rho, ntot = n_tot, nsim = mc$e2e) {
  d <- designs[[design]]
  r <- run_trial(pC, all_configs[[cfg]], n1, n_tot = ntot, rho = rho, gate = d$gate, w = d$w,
                 eta = get_eta(design, n1, rho), nsim = nsim, alpha = alpha, eta_pe = d$eta_pe)
  data.frame(cfg = cfg, design = design, n1 = n1, ntot = ntot, rho = rho, nsim = nsim,
             proceed = r$proceed, pcs = r$pcs, power = r$power, claim_opt = r$claim_opt,
             false_claim = r$false_claim, EN = r$EN)
}
cell_nsim <- function(cfg) if (cfg %in% names(null_configs)) mc$e2e_null else mc$e2e

n1_D <- c(30, n1_frontier)
ntot_D <- sort(unique(rd("e2e_runD.csv")$ntot))
matched0 <- rd("e2e_runD_matched.csv")

## Boundaries, computed serially so that the parallel workers share them
say("calibrating futility boundaries")
need <- unique(rbind(
  expand.grid(design = new_designs, n1 = 100, rho = rho_latent, stringsAsFactors = FALSE),
  expand.grid(design = new_designs, n1 = n1_D, rho = 0.3, stringsAsFactors = FALSE),
  data.frame(design = matched0$design, n1 = matched0$n1, rho = 0.3, stringsAsFactors = FALSE)))
for (i in seq_len(nrow(need))) get_eta(need$design[i], need$n1[i], need$rho[i])
eta_tab <- data.frame(key = ls(eta_cache), eta = unlist(mget(ls(eta_cache), envir = eta_cache)))
write_out(eta_tab, "e2e_pe_eta.csv")

## Block A
say("block A: n1 = 100")
cellsA <- merge(expand.grid(cfg = names(all_configs), design = new_designs, stringsAsFactors = FALSE),
                data.frame(rho = rho_latent))
runA <- do.call(rbind, papply(seq_len(nrow(cellsA)), function(i) {
  g <- cellsA[i, ]
  e2e_cell(g$cfg, g$design, 100, g$rho, nsim = cell_nsim(g$cfg))
}))
write_out(runA, "e2e_pe_runA.csv")

## Block D
say("block D: (n1, n_tot) grid")
cellsD <- rbind(
  expand.grid(cfg = "S1", design = new_designs, n1 = n1_D, ntot = ntot_D, stringsAsFactors = FALSE),
  expand.grid(cfg = "S3", design = "Cj_conc_pe", n1 = n1_D, ntot = ntot_D, stringsAsFactors = FALSE))
runD <- do.call(rbind, papply(seq_len(nrow(cellsD)), function(i) {
  g <- cellsD[i, ]
  e2e_cell(g$cfg, g$design, g$n1, 0.3, ntot = g$ntot)
}))
write_out(runD, "e2e_pe_runD.csv")

runD0 <- rd("e2e_runD.csv")
power_targets <- sort(unique(matched0$target))
matched <- do.call(rbind, lapply(unique(cellsD$cfg), function(s) {
  ds <- unique(cellsD$design[cellsD$cfg == s])
  do.call(rbind, lapply(power_targets, function(tp) {
    m <- min_EN(runD0[runD0$cfg == s & runD0$design == "M", ], tp)
    do.call(rbind, lapply(ds, function(g) {
      x <- min_EN(runD[runD$cfg == s & runD$design == g, ], tp)
      ratio <- x[["EN"]] / m[["EN"]]
      data.frame(cfg = s, target = tp, design = g,
                 n1 = round(x[["n1"]]), ntot = x[["ntot"]],
                 EN = round(x[["EN"]]), se_EN = round(x[["se_EN"]]),
                 at_bound = x[["bound"]], ratio = round(ratio, 3),
                 se_ratio = round(ratio * sqrt((x[["se_EN"]] / x[["EN"]])^2 +
                                               (m[["se_EN"]] / m[["EN"]])^2), 3))
    }))
  }))
}))
write_out(matched, "e2e_pe_runD_matched.csv")

## Block N: E[N] under N0 and N1 at each design's optimum (n1, n_tot)
say("block N: nulls at the run-D optima")
opt <- rbind(matched0[, c("cfg", "target", "design", "n1", "ntot")],
             matched[, c("cfg", "target", "design", "n1", "ntot")])
opt <- opt[!is.na(opt$n1), ]
for (i in seq_len(nrow(opt))) get_eta(opt$design[i], opt$n1[i], 0.3)
cellsN <- merge(cbind(row = seq_len(nrow(opt)), opt), data.frame(null = c("N0", "N1")))
nulls <- do.call(rbind, papply(seq_len(nrow(cellsN)), function(i) {
  g <- cellsN[i, ]
  r <- e2e_cell(g$null, g$design, g$n1, 0.3, ntot = g$ntot, nsim = mc$nulls)
  data.frame(cfg = g$cfg, target = g$target, design = g$design, n1 = g$n1, ntot = g$ntot,
             null = g$null, nsim = mc$nulls, proceed = r$proceed, false_claim = r$false_claim,
             EN = r$EN)
}))
write_out(nulls, "e2e_runD_nulls.csv")

## Figure fig:frontier with the two new curves
say("redrawing results/figures/sim_e2e_fig.pdf")
runD_all <- rbind(runD0, runD)
fig_design <- list(S1 = c("M", "Cu", "Cj_conc", "Cu_pe", "Cj_conc_pe"),
                   S3 = c("M", "Cc_S3", "Cj_conc", "Cj_conc_pe"))
fig_label <- list(M = quote("co-primary gate"),
                  Cu = quote(paste("composite, uniform ", bold(w))),
                  Cc_S3 = quote(paste("composite, ", bold(w)^"*", "(S3)")),
                  Cj_conc = quote(paste("composite, ", bold(w)^"*", (Theta[conc]))),
                  Cu_pe = quote(paste("uniform ", bold(w), " + per-endpoint check")),
                  Cj_conc_pe = quote(paste(bold(w)^"*", (Theta[conc]), " + per-endpoint check")))
fig_col <- c(M = "black", Cu = "#D55E00", Cc_S3 = "#D55E00", Cj_conc = "#0072B2",
             Cu_pe = "#D55E00", Cj_conc_pe = "#0072B2")
fig_lty <- c(M = 1, Cu = 2, Cc_S3 = 2, Cj_conc = 1, Cu_pe = 3, Cj_conc_pe = 3)
tgrid <- seq(0.60, 0.92, by = 0.01)
pdf(file.path(fig_dir, "sim_e2e_fig.pdf"), width = 9, height = 4.2)
op <- par(mfrow = c(1, 2), mar = c(3.8, 3.8, 2, 0.5), mgp = c(2.4, 0.7, 0))
for (s in names(fig_design)) {
  ds <- fig_design[[s]]
  env <- sapply(ds, function(g) sapply(tgrid, function(tp)
    min_EN(runD_all[runD_all$cfg == s & runD_all$design == g, ], tp)[["EN"]]))
  plot(NA, xlim = range(tgrid), ylim = range(env, na.rm = TRUE),
       xlab = "Target probability of co-primary claim",
       ylab = "Smallest expected total sample size", main = config_labels[s])
  for (g in ds) lines(tgrid, env[, g], col = fig_col[g], lty = fig_lty[g], lwd = 1.5)
  legend("topleft", legend = as.expression(fig_label[ds]), col = fig_col[ds],
         lty = fig_lty[ds], lwd = 1.5, bty = "n", cex = 0.85)
}
par(op); invisible(dev.off())
say("done", if (quick) "(quick mode)" else "")
