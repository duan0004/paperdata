# PRL Decisive Evidence Figure

**PDF**: `results/T2_NG15yr/figures/prl_decisive_evidence_figure.pdf`
**PNG**: `results/T2_NG15yr/figures/prl_decisive_evidence_figure.png`
**Data JSON**: `results/T2_NG15yr/figures/prl_decisive_evidence_figure_data.json`

## Panel A: Calibration Ablation

| Step | baseline | lnB | err | residual vs published |
|---|---|---:|---:|---:|
| local-template/local-KDE | fixed-gamma SMBHB | `-1.590` | `0.154` | `-5.633` |
| official-template/local-KDE | PTArcade BHB prior | `+0.687` | `0.147` | `-3.356` |
| official-template/official-density seed42 | PTArcade BHB prior | `+4.225` | `0.173` | `+0.182` |
| official-template/official-density robust | PTArcade BHB prior | `+4.520` | `0.185` | `+0.477` |

## Panel B: Matched-Scale Ranking

| Group | Model | lnB | plotted err |
|---|---|---:|---:|
| cosmology | sigw_gauss_ptarcade | `+4.520` | `0.237` |
| cosmology | sigw_delta_ptarcade | `+3.930` | `0.246` |
| curved SMBHB control | ecc_supp_fixed_gamma | `+3.839` | `0.158` |
| curved SMBHB control | broken_pl_fixed_gamma | `+3.656` | `0.217` |
| cosmology | super_ptarcade | `+3.558` | `0.205` |
| curved SMBHB control | env_fixed_gamma | `+3.453` | `0.165` |
| curved SMBHB control | env_free_gamma | `+3.360` | `0.255` |

## Interpretation

The figure shows that the absolute-evidence recovery requires the matched
official template, official density, and matched BHB-prior baseline.  Panel A
plots only the matched PTArcade BHB-prior ablation; the fixed-gamma local
diagnostic is retained in the data table.  On that matched scale, multiple curved-SMBHB
controls occupy the same evidence tier as the leading stochastic new-physics
templates.
