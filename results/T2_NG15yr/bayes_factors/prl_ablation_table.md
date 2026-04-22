# PRL Evidence-Ablation Table

**Date**: 2026-04-21  
**Purpose**: isolate which ingredients close the absolute Bayes-factor gap.

Published reference values are from the NANOGrav 15-year new-physics search
as recorded in project notes: SIGW-Gaussian `ln(57)=+4.043`, SIGW-delta
`ln(44)=+3.784`, and cosmic superstrings `ln(46)=+3.829`.

## SIGW-Gaussian Ablation

| Template | Density / likelihood | Baseline | lnB | Residual vs published |
|---|---|---|---:|---:|
| local semi-analytic SIGW-Gaussian | local scipy KDE from `hd_30f_fs.core` | fixed-gamma SMBHB | `-1.590 +/- 0.154` | `-5.633` |
| official PTArcade SIGW-Gaussian | local scipy KDE from `hd_30f_fs.core` | PTArcade BHB prior | `+0.687 +/- 0.147` | `-3.356` |
| official PTArcade SIGW-Gaussian | official PTArcade `ceffyl` density | PTArcade BHB prior, seed-42 run | `+4.225 +/- 0.173` | `+0.182` |
| official PTArcade SIGW-Gaussian | official PTArcade `ceffyl` density | PTArcade BHB prior, 5-config mean | `+4.520 +/- 0.185` | `+0.477` |

## Official-Density Leading Models

| Model | Density / likelihood | Baseline | robust mean lnB | residual vs published | max abs residual in sweep |
|---|---|---|---:|---:|---:|
| SIGW-Gaussian | official PTArcade `ceffyl` density | PTArcade BHB prior | `+4.520` | `+0.477` | `0.694` |
| SIGW-delta | official PTArcade `ceffyl` density | PTArcade BHB prior | `+3.930` | `+0.146` | `0.387` |
| Cosmic superstrings (`super.py`) | official PTArcade `ceffyl` density | PTArcade BHB prior | `+3.558` | `-0.271` | `0.491` |

## Conclusion

The ablation supports the PRL narrative:

1. Local KDE products preserve qualitative ranking but do not reproduce
   absolute evidence.
2. Replacing only the template is insufficient.
3. Matched official template, official density, and matched BHB-prior baseline
   restore the leading official-density evidences to sub-nat agreement.
