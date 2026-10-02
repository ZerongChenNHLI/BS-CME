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
##   eta      futility-boundary calibration
mc <- list(pcs = 10000, calib = 2000, cross = 20000, e2e = 10000, e2e_null = 40000, eta = 10000)
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
joint_specs <- list(
  all        = list(set = theta_sets$all,        floor = 0),
  concordant = list(set = theta_sets$concordant, floor = 0),
  all_floor  = list(set = theta_sets$all,        floor = w_floor)
)
wjoint <- papply(joint_specs, function(sp) {
  cons <- if (sp$floor > 0) function(w) all(w >= sp$floor - 1e-9) else NULL
  calibrate_weights(configs[sp$set], pC, n1 = 100, rho = 0.3, nsim = mc$calib,
                    constraint = cons)
})
names(wjoint) <- names(joint_specs)
w_joint_table <- data.frame(
  set = names(joint_specs),
  scenarios = sapply(joint_specs, function(sp) paste(sp$set, collapse = "+")),
  floor = sapply(joint_specs, `[[`, "floor"),
  t(sapply(wjoint, `[[`, "w")),
  mean_pcs_calibration = sapply(wjoint, `[[`, "pcs"))
names(w_joint_table)[4:7] <- paste0("w", 1:4)
write_out(w_joint_table, "weights_joint.csv")

cross_w <- list(Cu = w_uniform, Cj_all = wjoint$all$w, Cj_conc = wjoint$concordant$w,
                Cj_floor = wjoint$all_floor$w, Cc_S3 = w_star("S3", 0.3),
                Cc_S4 = w_star("S4", 0.3))
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

## Delta-method standard error of the interpolated n1: SE(PCS) divided by the
## slope of the PCS curve on the interpolation segment.
n1_se <- function(n1, pcs, target, nsim) {
  i <- which(pcs >= target)[1]
  if (is.na(i) || i == 1) return(NA_real_)
  slope <- (pcs[i] - pcs[i - 1]) / (n1[i] - n1[i - 1])
  sqrt(target * (1 - target) / nsim) / slope
}

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
##   Cc_S3  composite, weights calibrated on S3 alone (best case for S3)
##   Cc_S4  composite, weights calibrated on S4 alone (best case for S4)
designs <- list(
  M     = list(gate = "conjunctive", w = NULL),
  Cu    = list(gate = "composite", w = w_uniform),
  Cj    = list(gate = "composite", w = wjoint$all$w),
  Cc_S3 = list(gate = "composite", w = function(rho) w_star("S3", 0.3)),
  Cc_S4 = list(gate = "composite", w = function(rho) w_star("S4", 0.3))
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
                     target = proceed_null) {
  d <- designs[[design]]
  r <- run_trial(pC, theta, n1, n_tot = n_tot, rho = rho, gate = d$gate,
                 w = design_w(d, rho), eta = get_eta(design, n1, rho, target),
                 nsim = nsim, alpha = alpha, closure = closure)
  data.frame(design = design, n1 = n1, rho = rho, nsim = nsim, proceed = r$proceed,
             pcs = r$pcs, power = r$power, claim_opt = r$claim_opt,
             false_claim = r$false_claim, EN = r$EN)
}
cell_nsim <- function(cfg) if (cfg %in% names(null_configs)) mc$e2e_null else mc$e2e

## Pre-compute every boundary serially so that parallel workers share them.
say("calibrating futility boundaries")
need_eta <- unique(rbind(
  expand.grid(design = names(designs), n1 = 100, rho = rho_latent, target = proceed_null,
              stringsAsFactors = FALSE),
  expand.grid(design = c("M", "Cu", "Cc_S3"), n1 = n1_frontier, rho = 0.3, target = proceed_null,
              stringsAsFactors = FALSE),
  expand.grid(design = c("M", "Cu", "Cc_S3"), n1 = 100, rho = 0.3, target = proceed_null_grid,
              stringsAsFactors = FALSE)
))
for (i in seq_len(nrow(need_eta)))
  get_eta(need_eta$design[i], need_eta$n1[i], need_eta$rho[i], need_eta$target[i])

## Run A: equal stage-1 cost (n1 = 100), alternatives and nulls
say("run A")
runA_cells <- rbind(
  expand.grid(cfg = c("S1", "S2"), design = c("M", "Cu", "Cj"), stringsAsFactors = FALSE),
  expand.grid(cfg = "S3", design = c("M", "Cu", "Cj", "Cc_S3"), stringsAsFactors = FALSE),
  expand.grid(cfg = "S4", design = c("M", "Cu", "Cj", "Cc_S4"), stringsAsFactors = FALSE),
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

## Run C: frontier in n1 at rho_latent = 0.3
say("run C")
runC_cells <- rbind(
  expand.grid(cfg = "S1", design = c("M", "Cu"), n1 = n1_frontier, stringsAsFactors = FALSE),
  expand.grid(cfg = "S3", design = c("M", "Cc_S3"), n1 = n1_frontier, stringsAsFactors = FALSE)
)
runC <- do.call(rbind, papply(seq_len(nrow(runC_cells)), function(i) {
  g <- runC_cells[i, ]
  cbind(cfg = g$cfg, e2e_cell(configs[[g$cfg]], g$design, g$n1, 0.3))
}))
write_out(runC, "e2e_runC.csv")

## Matched power: interpolate each frontier at fixed power targets (Table tab:e2eCm)
## Delta-method standard errors: SE(power) on the interpolation segment divided
## by the slope of power in n1 gives SE(n1); E[N] is linear in n1 on the segment.
power_targets <- c(0.60, 0.65, 0.70, 0.75)
interp_at <- function(d, target) {
  d <- d[order(d$n1), ]
  i <- which(d$power >= target)[1]
  if (is.na(i) || i == 1) return(c(n1 = NA, EN = NA, se_n1 = NA, se_EN = NA))
  f <- (target - d$power[i - 1]) / (d$power[i] - d$power[i - 1])
  slope_p <- (d$power[i] - d$power[i - 1]) / (d$n1[i] - d$n1[i - 1])
  slope_EN <- (d$EN[i] - d$EN[i - 1]) / (d$n1[i] - d$n1[i - 1])
  se_n1 <- sqrt(target * (1 - target) / d$nsim[i]) / slope_p
  c(n1 = d$n1[i - 1] + f * (d$n1[i] - d$n1[i - 1]),
    EN = d$EN[i - 1] + f * (d$EN[i] - d$EN[i - 1]),
    se_n1 = se_n1, se_EN = se_n1 * abs(slope_EN))
}
matched <- do.call(rbind, lapply(c("S1", "S3"), function(s) {
  dM <- runC[runC$cfg == s & runC$design == "M", ]
  dC <- runC[runC$cfg == s & runC$design != "M", ]
  do.call(rbind, lapply(power_targets, function(tp) {
    m <- interp_at(dM, tp); cc <- interp_at(dC, tp)
    ratio <- cc["EN"] / m["EN"]
    data.frame(cfg = s, target = tp,
               n1_M = round(m["n1"]), EN_M = round(m["EN"]),
               n1_C = round(cc["n1"]), EN_C = round(cc["EN"]),
               se_n1_M = round(m["se_n1"], 1), se_n1_C = round(cc["se_n1"], 1),
               se_EN_M = round(m["se_EN"]), se_EN_C = round(cc["se_EN"]),
               ratio = round(ratio, 2),
               se_ratio = round(ratio * sqrt((cc["se_EN"] / cc["EN"])^2 + (m["se_EN"] / m["EN"])^2), 3))
  }))
}))
write_out(matched, "e2e_runC_matched.csv")

pdf(file.path(fig_dir, "sim_e2e_fig.pdf"), width = 9, height = 4.2)
op <- par(mfrow = c(1, 2), mar = c(3.8, 3.8, 2, 0.5), mgp = c(2.4, 0.7, 0))
for (s in c("S1", "S3")) {
  d <- runC[runC$cfg == s, ]
  plot(NA, xlim = range(d$EN), ylim = range(d$power), xlab = "Expected total sample size",
       ylab = "Probability of co-primary claim", main = config_labels[s])
  for (g in unique(d$design)) {
    dd <- d[d$design == g, ]; dd <- dd[order(dd$n1), ]
    col <- if (g == "M") "black" else "#0072B2"
    lines(dd$EN, dd$power, type = "o", pch = if (g == "M") 16 else 15, col = col)
  }
  legend("bottomright", c("co-primary gate", "composite gate"), pch = c(16, 15),
         col = c("black", "#0072B2"), lty = 1, bty = "n")
}
par(op); invisible(dev.off())
say("wrote results/figures/sim_e2e_fig.pdf")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
say("done", if (quick) "(quick mode)" else "")
