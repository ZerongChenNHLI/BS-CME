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
  bscme_methods.R        method code (outline only)
  scenarios.R            simulation configurations (outline only)
  run_simulations.R      simulation study (outline only)
```

## Compiling the manuscript

The class loads its fonts with `fontspec`, so XeLaTeX is required. From
`manuscript/`:

```
latexmk -xelatex main
```

On Overleaf, upload the contents of `manuscript/` and set Menu → Compiler → XeLaTeX.
