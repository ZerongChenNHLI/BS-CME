## =============================================================================
## Simulation configurations (manuscript Table "tab:scen" and Section e2e).
## =============================================================================

## Control seroresponse rates for serogroups A, C, W, Y (Yang et al 2024).
pC <- c(0.425, 0.497, 0.448, 0.434)
K <- length(pC)

## Latent correlations and their product-moment equivalents (phi ~ 0, 0.19, 0.41).
rho_latent <- c(0, 0.3, 0.6)

## Dose-effect configurations theta_{a,k} = p_{a,k} - p_{C,k}; rows = doses.
mk <- function(...) matrix(c(...), ncol = 4, byrow = TRUE)

configs <- list(
  S1 = mk(0.05, 0.05, 0.05, 0.05,
          0.10, 0.10, 0.10, 0.10,
          0.15, 0.15, 0.15, 0.15),
  S2 = mk(0.04, 0.06, 0.03, 0.05,
          0.12, 0.17, 0.10, 0.15,
          0.10, 0.15, 0.09, 0.13),
  S3 = mk(0.00, 0.16, 0.16, 0.16,
          0.07, 0.16, 0.16, 0.16,
          0.14, 0.16, 0.16, 0.16),
  S4 = mk(0.06, 0.06, 0.06, 0.12,
          0.12, 0.12, 0.12, 0.08,
          0.18, 0.18, 0.18, -0.02)
)
config_labels <- c(S1 = "S1 parallel", S2 = "S2 heterogeneous",
                   S3 = "S3 concentrated", S4 = "S4 trade-off")

## Null configurations for the end-to-end study. In N0-N2 every dose that is
## not fully effective is null on a common endpoint; in N3 dose a is null on
## endpoint a only, so no common null endpoint exists and per-endpoint closure
## does not control the false-claim rate (Appendix B).
null_configs <- list(
  N0 = matrix(0, 3, 4),                                     # global null
  N1 = cbind(configs$S1[, 1:3], 0),                         # endpoint 4 null
  N2 = rbind(matrix(0, 2, 4), configs$S1[3, ]),             # doses 1-2 null
  N3 = mk(0.00, 0.15, 0.15, 0.15,                           # staggered null
          0.15, 0.00, 0.15, 0.15,
          0.15, 0.15, 0.00, 0.15),
  N4 = mk(0.00, 0.25, 0.25, 0.25,                           # staggered null, large effects
          0.25, 0.00, 0.25, 0.25,
          0.25, 0.25, 0.00, 0.25)
)
null_labels <- c(N0 = "N0 global null", N1 = "N1 endpoint-4 null",
                 N2 = "N2 doses 1-2 null", N3 = "N3 staggered null",
                 N4 = "N4 staggered null, large effects")

## Parallel configuration with M doses spaced by 0.05 (Table "tab:msweep").
parallel_config <- function(M, step = 0.05) {
  matrix(rep(step * seq_len(M), times = 4), M, 4)
}

## Calibrated weights reported in the manuscript (rho_latent = 0.3, n1 = 100).
## run_simulations.R re-derives them; these are used when CALIBRATE = FALSE.
w_uniform <- rep(1 / 4, 4)
w_star_paper <- list(S3 = c(1, 0, 0, 0),
                     S4 = c(0.20, 0.20, 0.10, 0.50))

## Scenario sets for joint calibration (eq. weights) and the floor constraint
## cal{W} = {w : w_k >= w_floor for every k}, so that no component is ignored.
theta_sets <- list(all = c("S1", "S2", "S3", "S4"), concordant = c("S1", "S2", "S3"))
w_floor <- 0.10

## Futility targets Pr(proceed | global null) for the eta sensitivity run.
proceed_null_grid <- c(0.10, 0.20, 0.30)

## Stage-1 grid (14 values, 25 to 800) and end-to-end settings.
n1_grid <- c(25, 40, 50, 60, 75, 100, 125, 150, 200, 250, 300, 400, 600, 800)
n1_frontier <- c(40, 60, 80, 100, 120, 150, 200)
n_tot <- 350
alpha <- 0.025
pcs_target <- 0.8
proceed_null <- 0.20
