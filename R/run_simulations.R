## =============================================================================
## BS-CME simulation study: reproduces the tables and figures of Sections
## "saving" and "e2e" of manuscript/body.tex.
##
## Usage (from the repository root):
##   Rscript R/run_simulations.R            # full run, Monte Carlo sizes as in the paper
##   Rscript R/run_simulations.R --quick    # smoke test with ~10% of the replicates
##
## Options (environment variables):
##   BSCME_CORES=4       run cells in parallel with parallel::mclapply (not on Windows)
##   BSCME_CALIBRATE=0   use the weights reported in the paper instead of re-deriving them
##
## Output: results/*.csv and results/figures/{sim_fig,sim_e2e_fig}.pdf
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
cores <- as.integer(Sys.getenv("BSCME_CORES", "1"))
calibrate <- Sys.getenv("BSCME_CALIBRATE", "1") != "0"

here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else "R"
})
source(file.path(here, "bscme_methods.R"))
source(file.path(here, "scenarios.R"))

out_dir <- file.path(dirname(here), "results")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## Monte Carlo sizes (paper values; --quick divides them by ten)
##   pcs      selection-stage PCS curves          cross     cross-evaluation of weights
##   calib    weight calibration                  e2e       end-to-end cells, alternatives
##   e2e_null end-to-end cells, null configurations (error rates)
##   eta      futility-boundary calibration     probe    finite-sample level probe
##   wald     size of the one-sided Wald test
mc <- list(pcs = 10000, calib = 2000, cross = 20000, e2e = 10000, e2e_null = 40000, eta = 10000,
           probe = 200000, wald = 2000000)
if (quick) mc <- lapply(mc, function(x) x %/% 10)

## Reproducible parallelism: job i always runs on the i-th L'Ecuyer stream after
## the current seed, so the output is identical for every value of BSCME_CORES.
RNGkind("L'Ecuyer-CMRG")
set.seed(20260930)
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

## ---------------------------------------------------------------------------
## 0. Second-moment bracket (eq. bracket) and implied phi
## ---------------------------------------------------------------------------
phi <- sapply(rho_latent, function(r) phi_from_latent(mean(pC), r))
bracket <- t(sapply(phi, function(f) n1_ratio_bracket(K, f)))
write_out(data.frame(rho_latent, phi = round(phi, 3), round(bracket, 3),
                     V_K = round(min_normal_var(K), 3)),
          "bracket.csv")

## ---------------------------------------------------------------------------
## 1. Calibrated weights (eq. weights): per configuration and rho, n1 = 100
## ---------------------------------------------------------------------------
cells <- expand.grid(cfg = names(configs), rho = rho_latent,
                     stringsAsFactors = FALSE)
if (calibrate) {
  say("calibrating weights")
  wcal <- papply(seq_len(nrow(cells)), function(i) {
    calibrate_weights(list(configs[[cells$cfg[i]]]), pC, n1 = 100,
                      rho = cells$rho[i], nsim = mc$calib)$w
  })
} else {
  wcal <- lapply(seq_len(nrow(cells)), function(i)
    if (cells$cfg[i] %in% names(w_star_paper)) w_star_paper[[cells$cfg[i]]]
    else w_uniform)
}
w_table <- cbind(cells, do.call(rbind, wcal))
names(w_table)[3:6] <- paste0("w", 1:4)
w_table$concordant_uniform <- sapply(cells$cfg, function(s)
  is_concordant(configs[[s]], w_uniform))
w_table$concordant_calibrated <- mapply(function(s, w) is_concordant(configs[[s]], w),
                                        cells$cfg, wcal)
write_out(w_table, "weights_calibrated.csv")
w_star <- function(cfg, rho) unlist(w_table[w_table$cfg == cfg & w_table$rho == rho, 3:6])

## ---------------------------------------------------------------------------
## 1b. Joint calibration over scenario sets, as eq. (weights) prescribes, and
##     cross-evaluation of every weight vector on every configuration
##     (Tables tab:wjoint and tab:cross; rho = 0.3, n1 = 100)
## ---------------------------------------------------------------------------
say("joint calibration")
joint_free <- papply(theta_sets, function(set)
  calibrate_weights(configs[set], pC, n1 = 100, rho = 0.3, nsim = mc$calib))
names(joint_free) <- names(theta_sets)
## The floored solution is read from the same calibration sample as the free
## one: the best grid point with every weight >= w_floor. When the free optimum
## already satisfies the floor, the two coincide.
floor_best <- function(res, floor) {
  a <- res$all
  feas <- apply(a[, 1:K], 1, function(w) all(w >= floor - 1e-9))
  i <- which(feas)[which.max(a$pcs[feas])]
  list(w = as.numeric(a[i, 1:K]), pcs = a$pcs[i])
}
joint_specs <- list(
  all              = list(set = theta_sets$all,        floor = 0),
  concordant       = list(set = theta_sets$concordant, floor = 0),
  all_floor        = list(set = theta_sets$all,        floor = w_floor),
  concordant_floor = list(set = theta_sets$concordant, floor = w_floor)
)
wjoint <- list(all = joint_free$all, concordant = joint_free$concordant,
               all_floor = floor_best(joint_free$all, w_floor),
               concordant_floor = floor_best(joint_free$concordant, w_floor))
w_joint_table <- data.frame(
  set = names(joint_specs),
  scenarios = sapply(joint_specs, function(sp) paste(sp$set, collapse = "+")),
  floor = sapply(joint_specs, `[[`, "floor"),
  t(sapply(wjoint, `[[`, "w")),
  mean_pcs_calibration = sapply(wjoint, `[[`, "pcs"))
names(w_joint_table)[4:7] <- paste0("w", 1:4)
write_out(w_joint_table, "weights_joint.csv")

cross_w <- list(Cu = w_uniform, Cj_all = wjoint$all$w, Cj_conc = wjoint$concordant$w,
                Cj_floor = wjoint$all_floor$w, Cj_cfloor = wjoint$concordant_floor$w,
                Cc_S3 = w_star("S3", 0.3), Cc_S4 = w_star("S4", 0.3))
cross <- papply(names(configs), function(s) {
  r <- pcs_selection(pC, configs[[s]], 100, 0.3, weights = cross_w, nsim = mc$cross)
  data.frame(cfg = s, gate = names(r$pcs), pcs = unname(r$pcs), regret = unname(r$regret))
})
write_out(do.call(rbind, cross), "pcs_cross.csv")

## ---------------------------------------------------------------------------
## 2. Selection-stage PCS against n1 (Table tab:main, Figure sim_fig)
## ---------------------------------------------------------------------------
say("selection-stage PCS")
grid <- expand.grid(cfg = names(configs), rho = rho_latent, n1 = n1_grid,
                    stringsAsFactors = FALSE)
pcs <- papply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  r <- pcs_selection(pC, configs[[g$cfg]], g$n1, g$rho,
                     weights = list(Cu = w_uniform, Cc = w_star(g$cfg, g$rho)),
                     nsim = mc$pcs)
  c(r$pcs, setNames(r$regret, paste0("regret_", names(r$regret))))
})
pcs_curves <- cbind(grid, do.call(rbind, pcs))
write_out(pcs_curves, "pcs_curves.csv")

tab_main <- do.call(rbind, lapply(split(pcs_curves, list(pcs_curves$cfg, pcs_curves$rho)),
  function(d) {
    d <- d[order(d$n1), ]
    n <- sapply(c("M", "Cu", "Cc"), function(g) n1_at_target(d$n1, d[[g]], pcs_target))
    se <- sapply(c("M", "Cu", "Cc"), function(g) n1_se(d$n1, d[[g]], pcs_target, mc$pcs))
    # Entries beyond the grid are reported as bounds, as in the manuscript.
    ratio <- function(num, den) {
      top <- max(d$n1)
      if (is.finite(num) && is.finite(den)) return(sprintf("%.2f", num / den))
      if (!is.finite(num) && !is.finite(den)) return(NA_character_)
      if (!is.finite(den)) sprintf("<%.2f", num / top) else sprintf(">%.2f", top / den)
    }
    data.frame(cfg = d$cfg[1], rho = d$rho[1],
               n1_M = round(n["M"]), n1_Cu = round(n["Cu"]), n1_Cc = round(n["Cc"]),
               se_M = round(se["M"], 1), se_Cu = round(se["Cu"], 1), se_Cc = round(se["Cc"], 1),
               ratio_u = ratio(n["Cu"], n["M"]), ratio_c = ratio(n["Cc"], n["M"]))
  }))
tab_main <- tab_main[order(tab_main$cfg, tab_main$rho), ]
write_out(tab_main, "table_main_n1.csv")

pdf(file.path(fig_dir, "sim_fig.pdf"), width = 10, height = 8.5)
op <- par(mfrow = c(4, 3), mar = c(3.5, 3.5, 2, 0.5), mgp = c(2.2, 0.7, 0))
for (s in names(configs)) for (r in rho_latent) {
  d <- pcs_curves[pcs_curves$cfg == s & pcs_curves$rho == r, ]
  d <- d[order(d$n1), ]
  matplot(d$n1, d[, c("M", "Cu", "Cc")], type = "o", log = "x", ylim = c(0, 1),
          pch = c(16, 17, 15), lty = 1, col = c("black", "#D55E00", "#0072B2"),
          xlab = expression(n[1] ~ "per arm"), ylab = "PCS",
          main = sprintf("%s, rho = %.1f", config_labels[s], r))
  abline(h = pcs_target, lty = 2, col = "grey50")
  if (s == "S1" && r == 0)
    legend("bottomright", c("co-primary gate", "composite, uniform w",
                            "composite, calibrated w"),
           pch = c(16, 17, 15), col = c("black", "#D55E00", "#0072B2"), bty = "n")
}
par(op); invisible(dev.off())
say("wrote results/figures/sim_fig.pdf")

## ---------------------------------------------------------------------------
## 3. Dependence on the number of doses (Table tab:msweep)
## ---------------------------------------------------------------------------
say("M sweep")
msweep <- do.call(rbind, papply(2:5, function(M) {
  p <- sapply(n1_grid, function(n) pcs_selection(pC, parallel_config(M), n, 0.3,
                                                 list(Cu = w_uniform), nsim = mc$pcs)$pcs)
  nM <- n1_at_target(n1_grid, p["M", ], pcs_target)
  nC <- n1_at_target(n1_grid, p["Cu", ], pcs_target)
  data.frame(M = M, n1_M = round(nM), n1_C = round(nC),
             se_M = round(n1_se(n1_grid, p["M", ], pcs_target, mc$pcs), 1),
             se_C = round(n1_se(n1_grid, p["Cu", ], pcs_target, mc$pcs), 1),
             ratio = round(nC / nM, 2), delta_N = round((M - 1) * (nM - nC)))
}))
write_out(msweep, "table_msweep.csv")

## ---------------------------------------------------------------------------
## 4. End-to-end operating characteristics (Section e2e)
## ---------------------------------------------------------------------------

## A design: gate + weights; eta is calibrated under the global null.
##   M      conjunctive (co-primary) gate
##   Cu     composite, uniform weights
##   Cj     composite, weights calibrated jointly over S1-S4 (rho = 0.3)
##   Cj_conc composite, weights calibrated jointly over the concordant set S1-S3
##          subject to the floor w_k >= w_floor: the vector the design rule of
##          Section rule selects for that scenario set
##   Cc_S3  composite, weights calibrated on S3 alone (best case for S3)
##   Cc_S4  composite, weights calibrated on S4 alone (best case for S4)
designs <- list(
  M       = list(gate = "conjunctive", w = NULL),
  Cu      = list(gate = "composite", w = w_uniform),
  Cj      = list(gate = "composite", w = wjoint$all$w),
  Cj_conc = list(gate = "composite", w = wjoint$concordant_floor$w),
  Cc_S3   = list(gate = "composite", w = function(rho) w_star("S3", 0.3)),
  Cc_S4   = list(gate = "composite", w = function(rho) w_star("S4", 0.3))
)
design_w <- function(d, rho) if (is.function(d$w)) d$w(rho) else d$w

eta_cache <- new.env()
get_eta <- function(design, n1, rho, target = proceed_null) {
  key <- paste(design, n1, rho, target)
  if (is.null(eta_cache[[key]])) {
    d <- designs[[design]]
    eta_cache[[key]] <- calibrate_eta(pC, M = 3, n1 = n1, rho = rho, gate = d$gate,
                                      w = design_w(d, rho), target = target,
                                      nsim = mc$eta)
  }
  eta_cache[[key]]
}

e2e_cell <- function(theta, design, n1, rho, closure = "full", nsim = mc$e2e,
                     target = proceed_null, ntot = n_tot) {
  d <- designs[[design]]
  r <- run_trial(pC, theta, n1, n_tot = ntot, rho = rho, gate = d$gate,
                 w = design_w(d, rho), eta = get_eta(design, n1, rho, target),
                 nsim = nsim, alpha = alpha, closure = closure)
  data.frame(design = design, n1 = n1, ntot = ntot, rho = rho, nsim = nsim, proceed = r$proceed,
             pcs = r$pcs, power = r$power, claim_opt = r$claim_opt,
             false_claim = r$false_claim, EN = r$EN)
}
cell_nsim <- function(cfg) if (cfg %in% names(null_configs)) mc$e2e_null else mc$e2e

## Pre-compute every boundary serially so that parallel workers share them.
say("calibrating futility boundaries")
## Run D grid: stage-1 size per arm and per-arm total on the selected dose and control.
n1_D <- c(30, n1_frontier)
ntot_D <- c(350, 400, 450, 500, 600, 700, 800, 1000)
need_eta <- unique(rbind(
  expand.grid(design = names(designs), n1 = 100, rho = rho_latent, target = proceed_null,
              stringsAsFactors = FALSE),
  expand.grid(design = c("M", "Cu", "Cc_S3", "Cj_conc"), n1 = n1_D, rho = 0.3,
              target = proceed_null, stringsAsFactors = FALSE),
  expand.grid(design = c("M", "Cu", "Cc_S3"), n1 = 100, rho = 0.3, target = proceed_null_grid,
              stringsAsFactors = FALSE)
))
for (i in seq_len(nrow(need_eta)))
  get_eta(need_eta$design[i], need_eta$n1[i], need_eta$rho[i], need_eta$target[i])

## Run A: equal stage-1 cost (n1 = 100), alternatives and nulls
say("run A")
runA_cells <- rbind(
  expand.grid(cfg = c("S1", "S2"), design = c("M", "Cu", "Cj", "Cj_conc"), stringsAsFactors = FALSE),
  expand.grid(cfg = "S3", design = c("M", "Cu", "Cj", "Cj_conc", "Cc_S3"), stringsAsFactors = FALSE),
  expand.grid(cfg = "S4", design = c("M", "Cu", "Cj", "Cj_conc", "Cc_S4"), stringsAsFactors = FALSE),
  expand.grid(cfg = names(null_configs), design = names(designs), stringsAsFactors = FALSE)
)
runA_cells <- merge(runA_cells, data.frame(rho = rho_latent))
all_configs <- c(configs, null_configs)
runA <- do.call(rbind, papply(seq_len(nrow(runA_cells)), function(i) {
  g <- runA_cells[i, ]
  cbind(cfg = g$cfg, e2e_cell(all_configs[[g$cfg]], g$design, 100, g$rho,
                              nsim = cell_nsim(g$cfg)))
}))
write_out(runA, "e2e_runA.csv")

## Closure comparison (Table tab:closure): the full closed test of the manuscript
## ("full": all M*K dose-endpoint pairs) against the dose-level alternative
## ("union": other doses enter through q_a = max_k p1_{a,k}) and against closing
## each endpoint separately ("endpoint", not valid), at n1 = 100.
say("closure comparison")
closure_cells <- merge(
  expand.grid(cfg = c("S1", "S3", "N1", "N3", "N4"), design = c("M", "Cu"),
              closure = c("full", "union", "endpoint"), stringsAsFactors = FALSE),
  data.frame(rho = rho_latent))
closure <- do.call(rbind, papply(seq_len(nrow(closure_cells)), function(i) {
  g <- closure_cells[i, ]
  cbind(cfg = g$cfg, closure = g$closure,
        e2e_cell(all_configs[[g$cfg]], g$design, 100, g$rho, closure = g$closure,
                 nsim = cell_nsim(g$cfg)))
}))
write_out(closure, "e2e_closure.csv")

## Futility-target sensitivity (Table tab:eta): Pr(proceed | global null) of
## 0.10, 0.20, 0.30 at n1 = 100, rho = 0.3.
say("eta sensitivity")
eta_cells <- expand.grid(cfg = c("S1", "S3", "N1"), target = proceed_null_grid,
                         design = c("M", "Cu", "Cc_S3"), stringsAsFactors = FALSE)
eta_cells <- eta_cells[(eta_cells$cfg == "S3") == (eta_cells$design == "Cc_S3") |
                         eta_cells$design == "M", ]
etasens <- do.call(rbind, papply(seq_len(nrow(eta_cells)), function(i) {
  g <- eta_cells[i, ]
  cbind(cfg = g$cfg, target = g$target, eta = get_eta(g$design, 100, 0.3, g$target),
        e2e_cell(all_configs[[g$cfg]], g$design, 100, 0.3, nsim = cell_nsim(g$cfg),
                 target = g$target))
}))
write_out(etasens, "e2e_eta.csv")

## Run B: equal selection accuracy, each design at its own n1 from tab_main
say("run B")
runB_cells <- do.call(rbind, lapply(rho_latent, function(r) {
  t1 <- tab_main[tab_main$cfg == "S1" & tab_main$rho == r, ]
  t3 <- tab_main[tab_main$cfg == "S3" & tab_main$rho == r, ]
  data.frame(cfg = c("S1", "S1", "S3", "S3"), rho = r,
             design = c("M", "Cu", "M", "Cc_S3"),
             n1 = c(t1$n1_M, t1$n1_Cu, t3$n1_M, t3$n1_Cc))
}))
runB_cells <- runB_cells[is.finite(runB_cells$n1) & runB_cells$n1 < n_tot, ]
## Boundaries for the run-B sizes, serially, so that workers share them.
for (i in seq_len(nrow(runB_cells))) get_eta(runB_cells$design[i], runB_cells$n1[i], runB_cells$rho[i])
runB <- do.call(rbind, papply(seq_len(nrow(runB_cells)), function(i) {
  g <- runB_cells[i, ]
  cbind(cfg = g$cfg, e2e_cell(configs[[g$cfg]], g$design, g$n1, g$rho))
}))
write_out(runB, "e2e_runB.csv")

## Run C: frontier in n1 at rho_latent = 0.3 and fixed n_tot: co-primary gate
## against the composite gate with the weights the design rule selects.
say("run C")
runC_cells <- expand.grid(cfg = c("S1", "S3"), design = c("M", "Cj_conc"), n1 = n1_frontier,
                          stringsAsFactors = FALSE)
runC <- do.call(rbind, papply(seq_len(nrow(runC_cells)), function(i) {
  g <- runC_cells[i, ]
  cbind(cfg = g$cfg, e2e_cell(configs[[g$cfg]], g$design, g$n1, 0.3))
}))
write_out(runC, "e2e_runC.csv")

## Linear interpolation of a frontier (power against n1 at fixed n_tot) at a
## target claim probability. Delta-method standard errors: SE(power) on the
## interpolation segment divided by the slope of power in n1 gives SE(n1), and
## E[N] is linear in n1 on the segment. bound = 1 flags a target already met at
## the smallest n1 of the grid, where the frontier value is only an upper bound.
interp_at <- function(d, target) {
  d <- d[order(d$n1), ]
  i <- which(d$power >= target)[1]
  if (is.na(i)) return(c(n1 = NA, EN = NA, se_n1 = NA, se_EN = NA, bound = NA))
  if (i == 1) return(c(n1 = d$n1[1], EN = d$EN[1], se_n1 = NA, se_EN = NA, bound = 1))
  f <- (target - d$power[i - 1]) / (d$power[i] - d$power[i - 1])
  slope_p <- (d$power[i] - d$power[i - 1]) / (d$n1[i] - d$n1[i - 1])
  slope_EN <- (d$EN[i] - d$EN[i - 1]) / (d$n1[i] - d$n1[i - 1])
  se_n1 <- sqrt(target * (1 - target) / d$nsim[i]) / slope_p
  c(n1 = d$n1[i - 1] + f * (d$n1[i] - d$n1[i - 1]),
    EN = d$EN[i - 1] + f * (d$EN[i] - d$EN[i - 1]),
    se_n1 = se_n1, se_EN = se_n1 * abs(slope_EN), bound = 0)
}

## Run D: n_tot as a design parameter (Table tab:e2eD, Figure fig:frontier).
## Every (n1, n_tot) pair of the grid is simulated at rho_latent = 0.3 for the
## co-primary gate and for the composite gate with uniform weights (S1), the
## configuration-specific w*(S3) (S3) and the rule-selected weights (both);
## each design is then optimised over the pair for a target claim probability.
say("run D")
runD_cells <- rbind(
  expand.grid(cfg = "S1", design = c("M", "Cu", "Cj_conc"), n1 = n1_D, ntot = ntot_D,
              stringsAsFactors = FALSE),
  expand.grid(cfg = "S3", design = c("M", "Cc_S3", "Cj_conc"), n1 = n1_D, ntot = ntot_D,
              stringsAsFactors = FALSE))
runD <- do.call(rbind, papply(seq_len(nrow(runD_cells)), function(i) {
  g <- runD_cells[i, ]
  cbind(cfg = g$cfg, e2e_cell(configs[[g$cfg]], g$design, g$n1, 0.3, ntot = g$ntot))
}))
write_out(runD, "e2e_runD.csv")

## Smallest E[N] on the grid at a target: interpolate along n1 at each n_tot,
## then take the n_tot with the smallest interpolated E[N].
min_EN <- function(d, target) {
  per <- t(sapply(sort(unique(d$ntot)), function(nt)
    c(ntot = nt, interp_at(d[d$ntot == nt, ], target))))
  ok <- which(is.finite(per[, "EN"]))
  if (!length(ok)) return(c(ntot = NA, n1 = NA, EN = NA, se_n1 = NA, se_EN = NA, bound = NA))
  per[ok[which.min(per[ok, "EN"])], ]
}
power_targets <- c(0.80, 0.85, 0.90)
matchedD <- do.call(rbind, lapply(c("S1", "S3"), function(s) {
  ds <- unique(runD_cells$design[runD_cells$cfg == s])
  do.call(rbind, lapply(power_targets, function(tp) {
    m <- min_EN(runD[runD$cfg == s & runD$design == "M", ], tp)
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
write_out(matchedD, "e2e_runD_matched.csv")

## Finite-sample level probe (Table tab:level): the least favourable
## configurations with the futility rule switched off (eta = 0, every trial
## proceeds) and the co-primary gate, mc$probe trials per cell in chunks, so that
## the finite-sample size of the Wald test is visible above the bound of
## Theorem 1. "N4+" is N4 with the three non-null effects raised to 0.35.
say("level probe")
stag_null <- function(eff) { th <- matrix(eff, 3, K); th[cbind(1:3, 1:3)] <- 0; th }
stopifnot(isTRUE(all.equal(stag_null(0.25), null_configs$N4, check.attributes = FALSE)))
probe_cells <- data.frame(cfg = c("N4", "N4+", "N4+", "N4+"), eff = c(0.25, 0.35, 0.35, 0.35),
                          n1 = c(100, 100, 200, 40), ntot = c(350, 350, 350, 120),
                          rho = c(0.3, 0, 0, 0), stringsAsFactors = FALSE)
n_chunk <- 5
probe_jobs <- merge(cbind(cell = seq_len(nrow(probe_cells)), probe_cells),
                    data.frame(chunk = seq_len(n_chunk)))
probe_raw <- do.call(rbind, papply(seq_len(nrow(probe_jobs)), function(i) {
  g <- probe_jobs[i, ]
  r <- run_trial(pC, stag_null(g$eff), g$n1, n_tot = g$ntot, rho = g$rho,
                 gate = "conjunctive", eta = 0, nsim = mc$probe %/% n_chunk, alpha = alpha)
  data.frame(cell = g$cell, false_claim = r$false_claim, proceed = r$proceed)
}))
probe <- cbind(probe_cells, design = "M", nsim = (mc$probe %/% n_chunk) * n_chunk,
               proceed = tapply(probe_raw$proceed, probe_raw$cell, mean),
               false_claim = tapply(probe_raw$false_claim, probe_raw$cell, mean))
probe$se <- sqrt(probe$false_claim * (1 - probe$false_claim) / probe$nsim)
write_out(probe, "level_probe.csv")

## Size of the one-sided Wald test at nominal alpha for two arms of n subjects
## with a common success probability p (mc$wald pairs of binomial counts).
say("Wald size")
wald_cells <- expand.grid(n = c(40, 100, 150, 250), p = c(0.43, 0.50))
wald_size <- do.call(rbind, papply(seq_len(nrow(wald_cells)), function(i) {
  g <- wald_cells[i, ]
  xa <- rbinom(mc$wald, g$n, g$p); xc <- rbinom(mc$wald, g$n, g$p)
  data.frame(n = g$n, p = g$p, nsim = mc$wald,
             size = mean(wald_p(xa, xc, g$n, g$n) <= alpha))
}))
write_out(wald_size, "wald_size.csv")

## Figure fig:frontier: smallest E[N] on the (n1, n_tot) grid against the
## target claim probability, by design.
fig_label <- list(M = quote("co-primary gate"),
                  Cu = quote(paste("composite, uniform ", bold(w))),
                  Cc_S3 = quote(paste("composite, ", bold(w)^"*", "(S3)")),
                  Cj_conc = quote(paste("composite, ", bold(w)^"*", (Theta[conc]))))
fig_col <- c(M = "black", Cu = "#D55E00", Cc_S3 = "#D55E00", Cj_conc = "#0072B2")
fig_lty <- c(M = 1, Cu = 2, Cc_S3 = 2, Cj_conc = 1)
tgrid <- seq(0.60, 0.92, by = 0.01)
pdf(file.path(fig_dir, "sim_e2e_fig.pdf"), width = 9, height = 4.2)
op <- par(mfrow = c(1, 2), mar = c(3.8, 3.8, 2, 0.5), mgp = c(2.4, 0.7, 0))
for (s in c("S1", "S3")) {
  ds <- unique(runD_cells$design[runD_cells$cfg == s])
  env <- sapply(ds, function(g) sapply(tgrid, function(tp)
    min_EN(runD[runD$cfg == s & runD$design == g, ], tp)[["EN"]]))
  plot(NA, xlim = range(tgrid), ylim = range(env, na.rm = TRUE),
       xlab = "Target probability of co-primary claim",
       ylab = "Smallest expected total sample size", main = config_labels[s])
  for (g in ds) lines(tgrid, env[, g], col = fig_col[g], lty = fig_lty[g], lwd = 1.5)
  legend("topleft", legend = as.expression(fig_label[ds]), col = fig_col[ds],
         lty = fig_lty[ds], lwd = 1.5, bty = "n")
}
par(op); invisible(dev.off())
say("wrote results/figures/sim_e2e_fig.pdf")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
say("done", if (quick) "(quick mode)" else "")
