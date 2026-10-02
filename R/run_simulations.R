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
mc <- list(pcs = 4000, calib = 2000, e2e = 10000, eta = 10000)
if (quick) mc <- lapply(mc, function(x) x %/% 10)

papply <- function(X, FUN) {
  if (cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(X, FUN, mc.cores = cores, mc.set.seed = TRUE)
  } else lapply(X, FUN)
}
RNGkind("L'Ecuyer-CMRG")
set.seed(20260930)

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
## 2. Selection-stage PCS against n1 (Table tab:main, Figure sim_fig)
## ---------------------------------------------------------------------------
say("selection-stage PCS")
grid <- expand.grid(cfg = names(configs), rho = rho_latent, n1 = n1_grid,
                    stringsAsFactors = FALSE)
pcs <- papply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  pcs_selection(pC, configs[[g$cfg]], g$n1, g$rho,
                weights = list(Cu = w_uniform, Cc = w_star(g$cfg, g$rho)),
                nsim = mc$pcs)
})
pcs_curves <- cbind(grid, do.call(rbind, pcs))
write_out(pcs_curves, "pcs_curves.csv")

tab_main <- do.call(rbind, lapply(split(pcs_curves, list(pcs_curves$cfg, pcs_curves$rho)),
  function(d) {
    d <- d[order(d$n1), ]
    n <- sapply(c("M", "Cu", "Cc"), function(g) n1_at_target(d$n1, d[[g]], pcs_target))
    # Entries beyond the grid are reported as bounds, as in the manuscript.
    ratio <- function(num, den) {
      top <- max(d$n1)
      if (is.finite(num) && is.finite(den)) return(sprintf("%.2f", num / den))
      if (!is.finite(num) && !is.finite(den)) return(NA_character_)
      if (!is.finite(den)) sprintf("<%.2f", num / top) else sprintf(">%.2f", top / den)
    }
    data.frame(cfg = d$cfg[1], rho = d$rho[1],
               n1_M = round(n["M"]), n1_Cu = round(n["Cu"]), n1_Cc = round(n["Cc"]),
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
                                                 list(Cu = w_uniform), nsim = mc$pcs))
  nM <- n1_at_target(n1_grid, p["M", ], pcs_target)
  nC <- n1_at_target(n1_grid, p["Cu", ], pcs_target)
  data.frame(M = M, n1_M = round(nM), n1_C = round(nC), ratio = round(nC / nM, 2),
             delta_N = round((M - 1) * (nM - nC)))
}))
write_out(msweep, "table_msweep.csv")

## ---------------------------------------------------------------------------
## 4. End-to-end operating characteristics (Section e2e)
## ---------------------------------------------------------------------------

## A design: gate + weights; eta is calibrated under the global null.
designs <- list(
  M     = list(gate = "conjunctive", w = NULL),
  Cu    = list(gate = "composite", w = w_uniform),
  Cc_S3 = list(gate = "composite", w = function(rho) w_star("S3", 0.3)),
  Cc_S4 = list(gate = "composite", w = function(rho) w_star("S4", 0.3))
)
design_w <- function(d, rho) if (is.function(d$w)) d$w(rho) else d$w

eta_cache <- new.env()
get_eta <- function(design, n1, rho) {
  key <- paste(design, n1, rho)
  if (is.null(eta_cache[[key]])) {
    d <- designs[[design]]
    eta_cache[[key]] <- calibrate_eta(pC, M = 3, n1 = n1, rho = rho, gate = d$gate,
                                      w = design_w(d, rho), target = proceed_null,
                                      nsim = mc$eta)
  }
  eta_cache[[key]]
}

e2e_cell <- function(theta, design, n1, rho) {
  d <- designs[[design]]
  r <- run_trial(pC, theta, n1, n_tot = n_tot, rho = rho, gate = d$gate,
                 w = design_w(d, rho), eta = get_eta(design, n1, rho),
                 nsim = mc$e2e, alpha = alpha)
  data.frame(design = design, n1 = n1, rho = rho, proceed = r$proceed,
             pcs = r$pcs, power = r$power, false_claim = r$false_claim, EN = r$EN)
}

## Pre-compute every boundary serially so that parallel workers share them.
say("calibrating futility boundaries")
need_eta <- unique(rbind(
  expand.grid(design = names(designs), n1 = 100, rho = rho_latent, stringsAsFactors = FALSE),
  expand.grid(design = c("M", "Cc_S3"), n1 = n1_frontier, rho = 0.3, stringsAsFactors = FALSE),
  expand.grid(design = "Cu", n1 = n1_frontier, rho = 0.3, stringsAsFactors = FALSE)
))
for (i in seq_len(nrow(need_eta))) get_eta(need_eta$design[i], need_eta$n1[i], need_eta$rho[i])

## Run A: equal stage-1 cost (n1 = 100), alternatives and nulls
say("run A")
runA_cells <- rbind(
  expand.grid(cfg = c("S1", "S2"), design = c("M", "Cu"), stringsAsFactors = FALSE),
  expand.grid(cfg = "S3", design = c("M", "Cu", "Cc_S3"), stringsAsFactors = FALSE),
  expand.grid(cfg = "S4", design = c("M", "Cu", "Cc_S4"), stringsAsFactors = FALSE),
  expand.grid(cfg = names(null_configs), design = names(designs), stringsAsFactors = FALSE)
)
runA_cells <- merge(runA_cells, data.frame(rho = rho_latent))
all_configs <- c(configs, null_configs)
runA <- do.call(rbind, papply(seq_len(nrow(runA_cells)), function(i) {
  g <- runA_cells[i, ]
  cbind(cfg = g$cfg, e2e_cell(all_configs[[g$cfg]], g$design, 100, g$rho))
}))
write_out(runA, "e2e_runA.csv")

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
interp_at <- function(d, target) {
  d <- d[order(d$n1), ]
  i <- which(d$power >= target)[1]
  if (is.na(i) || i == 1) return(c(n1 = NA, EN = NA))
  f <- (target - d$power[i - 1]) / (d$power[i] - d$power[i - 1])
  c(n1 = d$n1[i - 1] + f * (d$n1[i] - d$n1[i - 1]),
    EN = d$EN[i - 1] + f * (d$EN[i] - d$EN[i - 1]))
}
matched <- do.call(rbind, lapply(c("S1", "S3"), function(s) {
  dM <- runC[runC$cfg == s & runC$design == "M", ]
  dC <- runC[runC$cfg == s & runC$design != "M", ]
  do.call(rbind, lapply(c(0.65, 0.70, 0.75, 0.80), function(tp) {
    m <- interp_at(dM, tp); cc <- interp_at(dC, tp)
    data.frame(cfg = s, target = tp, n1_M = round(m["n1"]), EN_M = round(m["EN"]),
               n1_C = round(cc["n1"]), EN_C = round(cc["EN"]),
               ratio = round(cc["EN"] / m["EN"], 2))
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
