# Codes_NGMMR_HRR

Replication code for the comment **"How Likely Is an Inflation Disaster? A Comment"** by Sophie Guilloux-Nefussi, Magali Marx, Sarah Mouabbi, and Jean-Paul Renne.

The comment revisits the mapping from option-implied risk-neutral inflation distributions to physical inflation tail probabilities in Hilscher, Raviv, and Reis (2025, hereafter HRR). The code implements two exercises:

1. Applying a model-consistent Epstein-Zin change of measure to the nominal risk-neutral transition matrices estimated by HRR.
2. Refitting parsimonious nominal risk-neutral transition matrices to HRR's option-implied probability targets, and then applying the same model-consistent change of measure.

The package also produces the three-date numerical example and the diagnostics on 5y5y probability targets used in the comment.

## How To Run

Open `Codes_NGMMR_HRR.Rproj` in RStudio and run:

```r
source("main.R")
```

By default, `main.R` uses the stored US and euro-area nominal-Q estimates and regenerates the figures and tables used in the comment. This is the fastest mode and typically takes about one minute on a standard computer.

To rerun the date-by-date estimation before producing the outputs, set:

```r
REESTIMATE <- TRUE
```

at the top of `main.R`. The estimation is parallelized with the `parallel` package; the number of cores is controlled by:

```r
N_CORES <- 8
```

With `REESTIMATE <- TRUE`, `main.R` reports progress by batches of dates and writes temporary `diagnostics_in_progress.csv` files in the relevant output folders.

## Software Requirements

The code is written in R and uses:

```r
library(R.matlab)
library(haven)
library(optimx)
library(parallel)
```

The remaining functions use base R.

## Data Sources

The package is self-contained. The folder `input/` contains upstream data objects from the HRR replication materials and public data files accompanying Hilscher, Raviv, and Reis (2025). The published HRR paper is available at <https://doi.org/10.1093/rfs/hhaf058>. These files are used as inputs and are not re-estimated here.

The HRR public data page is:

<https://r2rsquaredlse.github.io/web-inflationdisasters/>

The input files are:

- `input/results_est_101_US_monthquart.mat`
- `input/results_est_101_EZ_monthquart.mat`
- `input/dists_US.mat`
- `input/dists_EZ.mat`
- `input/USwestimates.dta`
- `input/EZwestimates.dta`
- `input/US_55tails_monthly.dta`
- `input/EZ_55tails_monthly.dta`

In the scripts:

- `results_est_101_*_monthquart.mat` contains HRR's estimated transition-matrix objects and model-implied probability targets.
- `dists_*.mat` contains HRR option-implied nominal and real risk-neutral distribution objects.
- `*westimates.dta` contains HRR reported physical probability series.
- `*_55tails_monthly.dta` contains HRR 5y5y tail-proxy inputs used for diagnostics and auxiliary targets.

No external data are downloaded by the scripts. The package reads the files in `input/`, regenerates derived estimates and figures, and writes outputs under `outputs/`.

## Main Outputs

Running `source("main.R")` regenerates:

- `outputs/comment_figures/figure_ratios.pdf`
- `outputs/comment_figures/HRR_model_consistent_P_gt4.pdf`
- `outputs/comment_figures/nominal_Q_tail_fit_HRR_NGMMR_gt4.pdf`
- `outputs/comment_figures/nominal_Q_tail_rmse_table.tex`
- `outputs/comment_figures/figure_IES_6plots_HRR_gt4.pdf`
- `outputs/comment_figures/figure_ratios_NGMMR.pdf`
- `outputs/hrr_nominal_real_q_diagnostics/5y5y_tail_bounds_from_5y_10y.pdf`
- `outputs/hrr_nominal_real_q_diagnostics/5y5y_gaussian_student_proxy_df5.pdf`
- `outputs/three_date_ez_example/three_date_ez_table_rows.tex`

When `REESTIMATE <- FALSE`, the code uses the stored estimates in:

- `outputs/US_nominal_Q_refit_smooth_row_poly2_finegrid_targets_hightailbin1p5_tailmoment2_moment5y5yW0p2_gauss5y5yRhodata_fwdBinW1/`
- `outputs/EZ_nominal_Q_refit_smooth_row_poly2_finegrid_targets_hightailbin1p5_tailmoment2_moment5y5yW0p2_gauss5y5yRhodata_fwdBinW1/`

Each folder contains:

- `all_results.rds`
- `ngmmr_Q_matrices.rds`
- `diagnostics.csv`

## Code Organization

All active scripts are stored in `scripts/`.

- `three_date_ez_example.R`: numerical three-date Epstein-Zin example and LaTeX table rows.
- `hrr_model_bins.R`: HRR transition-matrix construction and model-implied bin probabilities.
- `qn_parameterization.R`: active six-parameter nominal risk-neutral transition-matrix parameterization.
- `estimation_helpers.R`: optimization and date-loop helpers.
- `estimate_nominal_q.R`: date-by-date estimation of nominal-Q matrices for `AREA = "US"` or `AREA = "EZ"`.
- `pq_mapping_helpers.R`: nominal-to-real and real-to-physical probability mappings.
- `make_ratio_and_probability_figures.R`: ratio and physical-probability figures.
- `make_fit_diagnostics.R`: HRR-consistent-probability figures, nominal-Q fit diagnostics, and RMSE table.
- `make_5y5y_bounds_diagnostic.R`: Frechet-Hoeffding-bound diagnostic for the 5y5y high-inflation probability.
- `make_5y5y_student_proxy_diagnostic.R`: Gaussian and Student-t auxiliary 5y5y proxy diagnostic.

The active nominal-Q specification is the six-parameter `smooth_row_poly2` parameterization. It preserves HRR's zero pattern and uses quadratic destination-state logits whose slope and curvature vary smoothly with the current inflation state.

## Reproducibility Notes

Run the code from the project root, i.e., the folder containing `main.R`.

The default `REESTIMATE <- FALSE` mode is intended for fast reproduction of the figures and tables. The `REESTIMATE <- TRUE` mode documents and reruns the estimation step that produced the stored nominal-Q matrices.

## License

The original code and documentation in this repository are released under the MIT License; see `LICENSE`.

The files in `input/` are upstream HRR data objects included for replication. They are not covered by the MIT License for our code and remain subject to the terms of their original sources.
