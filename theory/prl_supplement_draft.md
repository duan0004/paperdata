# Supplemental Material: Calibrated PTA Evidence Favors Curved Spectra but Not a Unique Nanohertz Source

**Date**: 2026-04-23  
**Main draft**: `theory/paper_prl_submission.tex`  
**Status**: Markdown mirror of the compiled supplement source
`theory/prl_supplement.tex`.

This supplement records the data products, density convention, evidence
ablation, prior choices, curved-SMBHB control family, cross-PTA bridge
stress tests, convergence diagnostics, plateau diagnostics, and QMC/TI
cross-checks used by the PRL Letter.  The compiled submission source is
`theory/prl_supplement.tex`.

## S1. Data Products and Density Convention

The analysis uses the public NANOGrav 15-year data set (arXiv:2306.16213), the
NANOGrav 15-year noise analysis (arXiv:2306.16214), the NANOGrav new-physics
comparison paper (arXiv:2306.16219), official NANOGrav/PTArcade model files
from Zenodo DOI `10.5281/zenodo.8084351`, and the official PTArcade `ceffyl`
density product from Zenodo DOI `10.5281/zenodo.10495907`.

PTArcade model files return `h^2 Omega_GW`.  The residual-density conversion
used for official templates is

```text
rho^2 = (H0/h)^2 h^2Omega_GW / (8 pi^4 f^5 Tspan).
```

This conversion was checked against `ptarcade.models_utils.omega2cross` to
relative precision below `5e-16`.

## S2. Evidence-Ablation Check

| Template/density | Baseline | lnB | Residual |
|---|---|---:|---:|
| local template + local KDE | fixed-gamma SMBHB | `-1.590 +/- 0.154` | `-5.633` |
| official template + local KDE | PTArcade BHB prior | `+0.687 +/- 0.147` | `-3.356` |
| official template + official density, seed 42 | PTArcade BHB prior | `+4.225 +/- 0.173` | `+0.182` |
| official template + official density, five-config mean | PTArcade BHB prior | `+4.520 +/- 0.185` | `+0.477` |

The first row uses a fixed-gamma SMBHB baseline; the remaining rows use the
matched PTArcade BHB-prior baseline.

## S3. Official-Density Robustness

| Model | mean lnB | published | mean residual | max residual |
|---|---:|---:|---:|---:|
| SIGW-Gaussian | `+4.520` | `ln(57)=+4.043` | `+0.477` | `0.694` |
| SIGW-delta | `+3.930` | `ln(44)=+3.784` | `+0.146` | `0.387` |
| Cosmic superstrings | `+3.558` | `ln(46)=+3.829` | `-0.271` | `0.491` |

The correct main-text wording is sub-nat reproduction, not that every run is
within `0.5 nat`.

## S4. Priors and Model-Selection Rule

| Model/class | Parameters and priors | Rationale | Status |
|---|---|---|---|
| PTArcade BHB baseline | Gaussian prior on `(log10A, gamma)`, `mu=(-15.615,4.707)`, `C11=0.2787`, `C22=0.1242`, `C12=-0.00264` | matched NG15 PTArcade SMBHB baseline | fixed |
| SIGW-Gaussian | `log10 f_peak U[-11,-5]`, width `U[0.1,3]`, `log10 A U[-3,1]` | official PTArcade stochastic template | fixed |
| SIGW-delta | `log10 f_peak U[-11,-5]`, `log10 A U[-3,1]` | official PTArcade stochastic template; upper-bound variants reported separately | fixed |
| Cosmic superstrings | `log10 Gmu U[-14,-6]`, `log10 P U[-3,0]` | official PTArcade superstring template | fixed |
| Environmental turnover | `log10 A U[-18,-11]`, `log10 f_bend U[-9.8,-7.0]`, `kappa U[0.5,8]`, `gamma=13/3`; variants change `gamma` or bend range | SMBHB environmental-hardening control | fixed |
| Broken power law | `log10 A U[-18,-11]`, `log10 f_bend U[-10.0,-6.8]`, `delta U[0,6]`, `kappa U[0.5,8]`, `gamma=13/3` | phenomenological low-frequency suppression surrogate | fixed |
| Eccentricity-inspired | `log10 A U[-18,-11]`, `log10 f_e U[-10.0,-6.8]`, `beta U[0,3]`, `gamma=13/3` | phenomenological eccentricity-inspired suppression surrogate | fixed |

All curved-SMBHB controls were specified before the official-density family
sweep, and all tested variants are reported.  The Letter reports the
best-performing tested representative of each curvature class, not a search
over unreported curves.

## S5. Curved SMBHB Control Family

The environmental-turnover row is the astrophysical control.  The broken-power
law and eccentricity-inspired rows are phenomenological low-frequency
suppression controls, not population-synthesis source claims.

| Model | class | mean lnB | sample std | family TI residual |
|---|---|---:|---:|---:|
| eccentricity-inspired | ecc. | `+3.839` | `0.043` | `+0.056` |
| broken power law | BPL | `+3.656` | `0.154` | `+0.063` |
| env. fixed gamma=13/3 | env. | `+3.453` | `0.062` | `+0.017` |
| env. free gamma | env. | `+3.360` | `0.202` |  |
| env. broad bend | env. | `+2.992` | `0.142` |  |
| env. low bend | env. | `+2.666` | `0.103` |  |

At least two distinct curvature parameterizations overlap the leading
cosmological evidence tier within one nat.

The final column reports the thermodynamic-integration residual from the
family cross-check relative to the corresponding five-configuration dynesty
mean; blank entries were not included in that family cross-check.

## S6. Official Stochastic-Model Ranking

Only official PTArcade files defining `spectrum(f, ...)` were evaluated.  PBH
and ULDM files exposing `signal(toas, ...)` were skipped by design.

| Model | dim | lnB vs BHB prior |
|---|---:|---:|
| `sigw_gauss` | 3 | `+4.616 +/- 0.174` |
| `sigw_delta` | 2 | `+4.098 +/- 0.171` |
| `sigw_box` | 3 | `+3.830 +/- 0.182` |
| `super` | 2 | `+3.739 +/- 0.170` |
| `meta_ls` | 2 | `+3.677 +/- 0.181` |
| `dw_sm` | 4 | `+3.625 +/- 0.183` |
| `meta_l` | 2 | `+3.014 +/- 0.186` |
| `pt_bubble` | 6 | `+2.906 +/- 0.188` |
| `igw` | 3 | `+2.881 +/- 0.188` |
| `pt_sound` | 6 | `+2.238 +/- 0.192` |
| `dw_ds` | 4 | `+0.474 +/- 0.192` |
| `stable_m` | 1 | `-1.394 +/- 0.159` |
| `stable_k` | 1 | `-1.477 +/- 0.158` |
| `stable_n` | 1 | `-1.572 +/- 0.158` |
| `stable_c` | 1 | `-1.608 +/- 0.157` |

## S7. SIGW-Delta Prior-Boundary Sensitivity

| Prior | lnB | Delta lnB | best-fit log10 f_peak |
|---|---:|---:|---:|
| `log10 f_peak <= -5` | `+4.098 +/- 0.171` | `+0.000` | `-5.110` |
| `log10 f_peak <= -4` | `+4.300 +/- 0.171` | `+0.203` | `-4.148` |
| `log10 f_peak <= -3` | `+4.524 +/- 0.169` | `+0.426` | `-4.107` |

## S8-S11. Diagnostics

The supplement also includes the CAR covariance null calibration, chain
convergence and failed-chain handling, Sobol-QMC thermodynamic-integration
cross-check, and cumulative low-frequency-bin scan.

The environmental-row QMC/TI residual appears twice with different values by
design: `+0.064` is the five-scramble main cross-check, while `+0.017` is the
three-scramble family cross-check.  The cumulative-bin scan uses the
representative single configuration also used in the ranking table, so the
14-bin values need not equal the five-configuration means in the Letter.

## S12-S15. Cross-PTA Bridge and Version Diagnostics

The compiled supplement now includes:

- the full 11-template bridge model list, including family membership,
  official-template status, official-density-adapter status, and whether the
  model is a published NANOGrav calibration target;
- the sequential anchored bridge ablation:
  `NG15-off`, `NG15-off+PPTA-local`, `NG15-off+EPTA-local`, and `hybrid3`;
- the explicit statement that `hybrid3` is an anchored source-ranking stress
  test rather than a fully calibrated multi-PTA Bayes-factor scale;
- tested-representative family mixtures, not complete source-class
  marginalizations;
- the effective low-frequency slope projection,
  `gamma_eff ~ 3.7--4.1`, for the three curved-SMBHB controls across PPTA,
  NANOGrav, and EPTA windows;
- PTArcade version `1.1.5`, `ceffyl` version `1.41.2`, and the official-density
  grid-boundary diagnostic showing zero checked lower/upper-boundary fractions
  for SIGW-Gaussian, SIGW-delta, cosmic superstrings, eccentricity-inspired
  curvature, broken-power-law curvature, and environmental turnover.

The sequential ablation is saved at
`results/prl_reference_bridge/sequential_bridge_ablation.md`; the plateau
diagnostic is saved at
`results/T2_NG15yr/bayes_factors/prl_ceffyl_plateau_diagnostic.md`.

## S16. Minimal Local Verification Commands

Run from the project root:

```bash
python3 code/gwb_templates.py
python3 -m py_compile code/prl_*.py
python3 code/prl_decisive_evidence_figure.py
python3 code/prl_reference_bridge_pipeline.py p2seq --profile production
python3 code/prl_evidence_ti_qmc_crosscheck.py
python3 code/prl_ceffyl_plateau_diagnostic.py
python3 code/prl_package_static_gate.py
```

The package static gate writes `results/T2_NG15yr/prl_package_static_gate.json`
and `results/T2_NG15yr/prl_package_static_gate.md`.
