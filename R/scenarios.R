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

## Null configurations for the end-to-end study.
null_configs <- list(
  N0 = matrix(0, 3, 4),                                     # global null
  N1 = cbind(configs$S1[, 1:3], 0),                         # endpoint 4 null
  N2 = rbind(matrix(0, 2, 4), configs$S1[3, ])              # doses 1-2 null
)
null_labels <- c(N0 = "N0 global null", N1 = "N1 endpoint-4 null",
                 N2 = "N2 doses 1-2 null")

## Parallel configuration with M doses spaced by 0.05 (Table "tab:msweep").
parallel_config <- function(M, step = 0.05) {
  matrix(rep(step * seq_len(M), times = 4), M, 4)
}

## Calibrated weights reported in the manuscript (rho_latent = 0.3, n1 = 100).
## run_simulations.R re-derives them; these are used when CALIBRATE = FALSE.
w_uniform <- rep(1 / 4, 4)
w_star_paper <- list(S3 = c(0.95, 0, 0.05, 0),
                     S4 = c(0.20, 0.25, 0.05, 0.50))

## Stage-1 grid (14 values, 25 to 800) and end-to-end settings.
n1_grid <- c(25, 40, 50, 60, 75, 100, 125, 150, 200, 250, 300, 400, 600, 800)
n1_frontier <- c(40, 60, 80, 100, 120, 150, 200)
n_tot <- 350
alpha <- 0.025
pcs_target <- 0.8
proceed_null <- 0.20
