# PRL H4/H5/H6 Experiment Results

**Date**: 2026-04-21 14:27:02
**JSON**: `results/T2_NG15yr/bayes_factors/prl_H4_H5_H6_experiments.json`
**Baseline**: `smbhb_ptarcade_bhb_prior`

## H4: Astrophysical Curved-SMBHB Competitor

| Model | lnB vs BHB prior | lnB err | Best-fit parameters |
|---|---:|---:|---|
| `smbhb_env_fixed_gamma` | `+3.814` | `0.179` | `{'log10_A': -14.138147212059216, 'log10_fbend': -7.813957443833466, 'kappa': 1.848396628934965}` |
| `smbhb_env_free_gamma` | `+3.507` | `0.181` | `{'log10_A': -14.437017302221175, 'log10_fbend': -7.7713805179855955, 'kappa': 3.9772655555278096, 'gamma': 6.523658602706696}` |

## H5: Official-Density Stochastic-Background Model Ranking

Only PTArcade model files defining `spectrum(f, ...)` were evaluated.
PBH and ULDM `signal(toas, ...)` models were skipped by design.

| Rank | Model | dim | lnB vs BHB prior | lnB err | published lnB | residual |
|---:|---|---:|---:|---:|---:|---:|
| 1 | `sigw_gauss` | 3 | `+4.616` | `0.174` | `+4.043` | `+0.573` |
| 2 | `sigw_delta` | 2 | `+4.098` | `0.171` | `+3.784` | `+0.313` |
| 3 | `sigw_box` | 3 | `+3.830` | `0.182` |  |  |
| 4 | `super` | 2 | `+3.739` | `0.170` | `+3.829` | `-0.090` |
| 5 | `meta_ls` | 2 | `+3.677` | `0.181` |  |  |
| 6 | `dw_sm` | 4 | `+3.625` | `0.183` |  |  |
| 7 | `meta_l` | 2 | `+3.014` | `0.186` |  |  |
| 8 | `pt_bubble` | 6 | `+2.906` | `0.188` | `+2.890` | `+0.016` |
| 9 | `igw` | 3 | `+2.881` | `0.188` |  |  |
| 10 | `pt_sound` | 6 | `+2.238` | `0.192` | `+1.308` | `+0.930` |
| 11 | `dw_ds` | 4 | `+0.474` | `0.192` |  |  |
| 12 | `stable_m` | 1 | `-1.394` | `0.159` |  |  |
| 13 | `stable_k` | 1 | `-1.477` | `0.158` |  |  |
| 14 | `stable_n` | 1 | `-1.572` | `0.158` |  |  |
| 15 | `stable_c` | 1 | `-1.608` | `0.157` |  |  |

## H6: SIGW-Delta Prior-Boundary Sensitivity

| Prior label | lnB vs BHB prior | lnB err | delta lnB vs original | Best-fit parameters |
|---|---:|---:|---:|---|
| `sigw_delta_prior_fmax_-3` | `+4.524` | `0.169` | `+0.426` | `{'log10_f_peak': -4.106950208272414, 'log10_A': 0.8228830175184667}` |
| `sigw_delta_prior_fmax_-4` | `+4.300` | `0.171` | `+0.203` | `{'log10_f_peak': -4.148210088842735, 'log10_A': 0.7880985220284034}` |
| `sigw_delta_prior_fmax_-5` | `+4.098` | `0.171` | `+0.000` | `{'log10_f_peak': -5.109932388995126, 'log10_A': -0.030234585539879788}` |

## Interpretation Placeholder

Use this table to decide whether the PRL claim remains merely a provenance
audit or becomes a broader source-identification systematic: if the curved
SMBHB model is competitive with the leading cosmological templates, the PRL
story should foreground curved-SMBHB source-identification degeneracy.
