# PRL Official-Density Robustness Sweep

**Date**: 2026-04-21  
**Script**: `code/prl_official_density_robustness.py`  
**JSON output**: `results/T2_NG15yr/bayes_factors/prl_official_density_robustness.json`  
**Log**: `logs/prl_official_density_robustness_20260421_135224.log`

## Setup

Evidence calculations use:

- official PTArcade `ceffyl` density product, Zenodo DOI `10.5281/zenodo.10495907`;
- official NANOGrav/PTArcade model files, Zenodo DOI `10.5281/zenodo.8084351`;
- 14 low-frequency bins;
- PTArcade BHB-prior SMBHB baseline;
- `dynesty` static nested sampling.

Sampler configurations:

| label | nlive | dlogz | seed |
|---|---:|---:|---:|
| nlive500_seed42 | 500 | 0.1 | 42 |
| nlive500_seed7 | 500 | 0.1 | 7 |
| nlive500_seed123 | 500 | 0.1 | 123 |
| nlive1000_seed42 | 1000 | 0.05 | 42 |
| nlive1000_seed7 | 1000 | 0.05 | 7 |

## Aggregate Result

Residuals are relative to the published NANOGrav new-physics values recorded
in the project notes: `ln(57)`, `ln(44)`, and `ln(46)`.

| Model | mean lnB | sample std | mean nested lnB err | mean residual | max abs residual |
|---|---:|---:|---:|---:|---:|
| SIGW-Gaussian | `+4.520` | `0.185` | `0.148` | `+0.477` | `0.694` |
| SIGW-delta | `+3.930` | `0.199` | `0.145` | `+0.146` | `0.387` |
| Cosmic superstrings (`super.py`) | `+3.558` | `0.145` | `0.145` | `-0.271` | `0.491` |

## Interpretation For PRL

The correct wording is **sub-nat reproduction**, not "every run is below
0.5 nat".  Across five seed/live-point configurations, all three leading
official-density models remain within `0.70 nat` of the published log-Bayes
factors, and the mean residuals are within `0.5 nat`.

This strengthens the PRL claim because the result is not a single-seed
coincidence.  The remaining `~0.15-0.20 nat` sample scatter is comparable to
the nested-sampling uncertainty and should be reported as part of the evidence
calibration budget.
