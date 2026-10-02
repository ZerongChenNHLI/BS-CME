# BS-CME
Bayesian seamless II/III trial with endpoint switch

A Bayesian seamless phase II/III dose-optimisation design that selects the dose
on a pre-specified weighted composite of $K$ binary endpoints and confirms the
selected dose on the full co-primary criterion, with stage-1 data included.
Manuscript in preparation for *Biometrical Journal*.

## Repository layout

```
manuscript/              LaTeX source (Wiley NJD v5 template)
  main.tex               preamble, title page, abstract
  body.tex               main text
  refs.bib               references
  figures/               sim_fig.pdf, sim_e2e_fig.pdf
  Fonts/<Family>/        fonts loaded by WileyNJDv5.cls via ./Fonts/<Family>/
  WileyNJDv5.cls, wileyNJD-AMA.bst, NJDnatbib.sty, LETTERSP.STY, latexmkrc
R/
  bscme_methods.R        method code (model, gates, final analysis, calibration)
  scenarios.R            simulation configurations S1-S4, N0-N2
  run_simulations.R      reproduces the simulation tables and figures
results/                 output of run_simulations.R (CSV tables, figures)
```

## Compiling the manuscript

The class loads its fonts with `fontspec`, so XeLaTeX is required. From
`manuscript/`:

```
latexmk -xelatex main
```

On Overleaf, upload the contents of `manuscript/` and set Menu → Compiler → XeLaTeX.

## Running the code

Base R only (≥ 4.0); no packages need installing.

```
Rscript R/run_simulations.R --quick           # smoke test, ~1 min
BSCME_CORES=4 Rscript R/run_simulations.R     # full Monte Carlo sizes, ~5 min on 4 cores
```

`BSCME_CALIBRATE=0` uses the calibrated weights reported in the paper instead of
re-deriving them. Outputs land in `results/`; copy `results/figures/*.pdf` to
`manuscript/figures/` to update the paper's figures.

| Output | Manuscript |
| --- | --- |
| `bracket.csv` | second-moment bracket, eq. (bracket) |
| `weights_calibrated.csv` | calibrated weights, eq. (weights) |
| `table_main_n1.csv`, `figures/sim_fig.pdf` | Table tab:main, Figure fig:pcs |
| `table_msweep.csv` | Table tab:msweep |
| `e2e_runA.csv` | Tables tab:e2eA and tab:e2eN |
| `e2e_runB.csv` | Table tab:e2eB |
| `e2e_runC.csv`, `e2e_runC_matched.csv`, `figures/sim_e2e_fig.pdf` | Tables tab:e2eC, tab:e2eCm, Figure fig:frontier |

### Using the method functions

```r
source("R/bscme_methods.R"); source("R/scenarios.R")

# Stage 1 of one trial: n1 = 100 per arm, configuration S1, latent rho = 0.3
x1  <- simulate_counts(config_cell_probs(pC, configs$S1, rho = 0.3), n = 100, nsim = 1)
mom <- posterior_moments(x1)                      # exact Dirichlet moments
gs  <- gate_stats(mom, w = w_uniform)             # composite and conjunctive z
eta <- calibrate_eta(pC, M = 3, n1 = 100, rho = 0.3, w = w_uniform)
select_dose(gs$zC, eta)                           # selected dose and go/no-go

# Operating characteristics of the full design
run_trial(pC, configs$S1, n1 = 100, n_tot = 350, rho = 0.3,
          gate = "composite", w = w_uniform, eta = eta, nsim = 2000)[1:5]
```
