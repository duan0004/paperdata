# Data Manifest

## Included Inputs

| Path | Purpose |
|---|---|
| `data/NG15yr/PTArcade_models_1.0.0/` | Official PTArcade model files used for the stochastic template evidence calculations. |
| `data/NG15yr/PTArcade_ceffyl_0.2.0/` | Official PTArcade `ceffyl` free-spectrum density product used for absolute Bayes factors. |
| `data/NG15yr/data_release/` | NANOGrav 15-year release metadata included with the local data checkout. |
| `data/NG15yr/tutorials/data/` | NANOGrav 15-year tutorial data used by the ENTERPRISE/PTMCMC validation scripts. |
| `data/NG15yr/tutorials/presampled_cores/` | Public pre-sampled cores used for free-spectrum posterior and covariance diagnostics. |

The nested `.git` directory from the original NANOGrav data checkout is not
included.  Large local MCMC chain directories under `results/T2_NG15yr/chains*`
are not included; summary JSON diagnostics are included instead.

## Included Derived Results

| Path | Purpose |
|---|---|
| `results/T2_NG15yr/bayes_factors/` | Evidence tables, prior-boundary scans, curved-SMBHB family results, QMC/TI cross-checks, and systematic envelope. |
| `results/T2_NG15yr/covariance/` | Covariance extraction products and CAR/null-calibration diagnostics. |
| `results/T2_NG15yr/figures/` | PRL main figure and cumulative-bin figure. |
| `results/T2_NG15yr/T2_*json` | Compact MCMC/posterior validation summaries. |

## Size Policy

The repository includes files needed for code/data review and compact
reproduction, while avoiding local chain dumps and transient partial files.
No included file is larger than the standard GitHub 100 MB file limit.
