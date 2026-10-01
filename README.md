# BAOS — Benchmark-Aware Objective Selection

Toolbox and reproduction code for:

> J.G.M. Mukuna, *When does a benchmark-relative objective change the optimum? Theory and a
> resolution-aware independence test for multi-objective design*, Structural and Multidisciplinary
> Optimization (submitted, 2026).

The toolbox decides, **before optimization**, whether a benchmark-relative objective `A = P/B`
(performance divided by a theoretical or attainable limit) adds a genuinely new optimization
dimension, and quantifies how far it moves the optimum and the Pareto set.

## Repository structure

| Folder | Content |
|---|---|
| `toolbox/` | Reusable BAOS toolbox: BOIT (global and local classes, confidence intervals), objective selection, exact lexicographic optima, optimum-shift measures, Pareto utilities, Theorem 1 check, worked example and regression tests |
| `paper_reproduction/pipeline/` | Complete pipeline of the paper (models, optimizers, surrogates, analyses), all computed results (`results/`, CSV and MAT) and the figure script |
| `paper_reproduction/revision_studies/` | Scripts and data of the revision studies: two-tier classification, noise and correlated-noise studies, graded-independence calibration, high-dimensional tests, Monte Carlo regrets, composite case |

## Requirements

GNU Octave 6 or later, or MATLAB R2016b or later. No toolboxes required. Figures are rendered
with Python 3 (NumPy, pandas, Matplotlib) from the saved CSV results.

All results were produced in GNU Octave 8.4. MATLAB compatibility is a design goal but has not
been tested; the reproduction pipeline uses a few Octave-only commands (`printf`, `fflush`,
`save('-v7', ...)`) that must be adapted for MATLAB. The `toolbox/` folder runs in both.

## Quick start (toolbox)

```matlab
run('toolbox/setup_baos.m')
test_boit                     % known-answer regression tests
example_custom_problem        % worked example on a toy heat engine
```

## Reproducing the paper

```matlab
cd paper_reproduction/pipeline
run_all                       % stages st01-st09; master seed 20260930
```
then `python3 plot_results_v2.py` for the figures. The revision studies are run with the scripts in
`paper_reproduction/revision_studies/` (see the comments at the top of each script).

## Citation and licence

Please cite the article above (see `CITATION.cff`). Released under the MIT licence.
