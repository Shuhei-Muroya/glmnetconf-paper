# Reproducibility materials for "Data-driven configuration tuning of glmnet for balancing accuracy and computational efficiency"

Shuhei Muroya and Kei Hirose

This repository contains the code, the generated summary dataset, the trained
network and the raw simulation outputs behind every figure and table of the
paper. The proposed method itself is implemented in the R package
[`glmnetconf`](https://github.com/Shuhei-Muroya/glmnetconf).

## Structure

Every directory corresponds to one part of the paper. Each one holds the
scripts, the stored outputs (`results/`) and the figures of the paper
(`figures/`).

| Directory | Paper | Scripts | Stored outputs |
|---|---|---|---|
| `R/` | shared functions | `covariance.R`, `data_generation.R`, `lasso_comparison.R`, `load_results.R` | |
| `01_summary_dataset/` | Sec. 3.2.1, App. B | `recipe.R` → `run_generation.R` → `build_summary_dataset.R` | `results/<part>/res_chunk_*.rds`, `summary_dataset.csv.gz` |
| `02_glmnet_mlp/` | Sec. 3.2.2, App. C | `train_glmnet_mlp.py`, `compare_activation.py`, `export_weights.py`, `acquisition.py` | `results/eig{10,5,1}/`, `results/eigenvalue_*.csv` |
| `03_illustrative_example/` | Fig. 1, Fig. 3 | `fig1_fig3.R` | `figures/fig1*.pdf`, `figures/fig3*.pdf` |
| `04_simulation/` | Sec. 4.1, Figs. 4–7 | `recipe.R` → `run_simulation.R` → `plot_figures.R` | `results/res_chunk_*.rds` |
| `05_compressed_sensing/` | Sec. 4.2, Fig. 8 | `run_reconstruction.R` → `plot_reconstruction.py` | `results/reconst_*.csv`, `results/rmse_scaled.csv` |
| `06_real_data/` | Sec. 4.3, Table 1 | `prepare_*.R` → `run_*.R` → `table1.R` | `data/*.csv`, `results/*.csv` |
| `A_configuration_effect/` | App. A, Figs. 9–10 | `recipe.R` → `run_simulation.R` → `plot_figures.R` | `results/res_chunk_*.rds` |
| `D_sensitivity/` | App. D, Figs. 12–15 | `recipe_T_hope.R`, `recipe_K.R` → `run_simulation.R` → `plot_T_hope.R`, `plot_K.R` | `results/{T_hope,K}/res_chunk_*.rds` |
| `E_lars_time_prediction/` | App. E.2, Fig. 17 | `run_forecast.R` → `plot_figure17.R` | `results/forecast_lars_all.csv` |

## How to run

All R scripts are run with the repository root as the working directory
(open `glmnetconf-paper.Rproj`, or `setwd()` to the root), e.g.

```r
source("04_simulation/plot_figures.R")
```

Python scripts in `02_glmnet_mlp/` are run from that directory; the one in
`05_compressed_sensing/` is run from the repository root.

* `run_*.R` / `run_generation.R` reproduce the raw outputs from scratch. They
  use `doParallel` with 20 cores and take hours to days; their outputs are
  already stored in `results/`, so this step can be skipped.
* `plot_*.R`, `build_summary_dataset.R` and `table1.R` only read the stored
  outputs and run in minutes.

## The summary dataset (Sec. 3.2.1)

`01_summary_dataset/summary_dataset.csv.gz` (gzip-compressed CSV, 25 MB;
231 MB uncompressed) has 567,000 rows
(63,000 datasets × 9 configurations) and the columns

```
n, p, sample_ev_top1..10, sample_ev_bottom1..10, nlambda, thresh, time_glmnet, SPE_glmnet
```

i.e. (N, p, γ₁..γ₁₀, γ₋₁..γ₋₁₀, n_λ, τ, T_glmnet, SPE). It is built from the
raw outputs in `01_summary_dataset/results/` by `build_summary_dataset.R`:

| `results/` | design distribution | SNR | datasets |
|---|---|---|---|
| `xdist1_normal/` | multivariate normal | 0.05–6 (5 values) | 11,250 |
| `xdist2_lognormal/` | log-normal | 0.05–6 | 11,250 |
| `xdist3_mixture/` | mixture of normals | 0.05–6 | 11,250 |
| `xdist4_nonlinear_factor/` | nonlinear factor model | 0.05–6 | 11,250 |
| `high_snr/` | all four | 100, 500 | 18,000 |

Each `res_chunk_*.rds` additionally keeps, for every dataset, the full
simulation setting (covariance structure, β structure, SNR, seed), the
eigenvalues of the true covariance matrix, and the results of glmnet
(default) and LARS.

## Trained network (Sec. 3.2.2, App. C)

`02_glmnet_mlp/results/eig10/` is the final glmnet-MLP (top/bottom 10
eigenvalues, hidden layers 49–63–32): the PyTorch checkpoint, the Optuna study
(`*.db`, 500 trials), the input/output scalers and the weights exported as
CSV (`model_params_glmnet_eig10/`). These CSV files and scalers are identical
to the ones shipped in `glmnetconf` (`inst/extdata/`). `eig5/` and `eig1/`
are the networks of the eigenvalue-feature comparison in Appendix C.

## Verification

The following checks were made on the stored outputs:

* `build_summary_dataset.R` rebuilds `summary_dataset.csv.gz` exactly (same rows
  in the same order) from `01_summary_dataset/results/`.
* The recipes reproduce the settings and seeds of all stored datasets in
  `01_summary_dataset`, `04_simulation`, `A_configuration_effect` and
  `D_sensitivity`.
* Regenerating datasets from their recipes and seeds: for `xdist1_normal`
  the stored eigenvalue features, true eigenvalues and noise levels were
  reproduced exactly for all sampled datasets. For the parts generated on the
  server, the true eigenvalues were always reproduced, and datasets with
  Σ = I were reproduced exactly on a different machine (macOS arm64,
  LAPACK 3.11.0), but most others were not: `MASS::mvrnorm()` uses an
  eigendecomposition of Σ whose eigenvectors depend on the BLAS/LAPACK build.
  Bit-identical regeneration therefore requires the environment of
  Appendix B.3 (x86_64, reference BLAS/LAPACK 3.12.0).
* The plotting scripts reproduce Figs. 4–7, 9, 10, 12–15 and 17
  pixel-identically from the stored outputs; `table1.R` reproduces Table 1;
  the RMSE values of Fig. 8 are recomputed from the stored reconstructions.

## Notes

* `xdist1_normal/` was generated first, with `gen_data_flexible_v1()`
  (no column scaling of X; noise variance βᵀΣβ / SNR) and a recipe built with
  the normal design only (`master_recipe_xdis1` in `recipe.R`).
  `run_generation.R` selects these automatically for this part.
* In the summary dataset, the configurations are
  τ ∈ {1e-7, 1e-8, 1e-9} and n_λ ∈ {⌊0.5p⌋, p, ⌊1.5p⌋} (`run_generation.R`).
* Computation times depend on the machine (Appendix B.3: AMD EPYC 7763,
  R 4.3.3, glmnet 4.1.8, lars 1.3, reference BLAS/LAPACK).

## Environment

R (4.3.x) packages: `glmnetconf`, `glmnet`, `lars`, `MASS`, `Matrix`,
`igraph`, `doParallel`, `foreach`, `dplyr`, `tidyr`, `purrr`, `ggplot2`,
`patchwork` (>= 1.2.0), `LimeSoDa`.
Python packages: see `02_glmnet_mlp/requirements.txt`.
