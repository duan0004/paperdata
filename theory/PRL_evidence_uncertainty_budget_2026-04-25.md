# PRL Evidence Uncertainty Budget

Date: 2026-04-25

## Scope

This note consolidates the order-unity evidence gaps and the calibration /
numerical / prior-boundary systematics relevant to the current PRL manuscript.
It is a support table for interpretation.  It is not a new model search and not
a replacement for the official-density evidence tables.

## Main Evidence Gaps

| quantity | size | source |
|---|---:|---|
| NG15 official-density gap: leading SIGW-Gaussian model minus best curved-SMBHB control | 0.681 nat | `theory/PRL_submission_package_manifest.md` |
| Hybrid3 tested-representative family gap: SIGW-like family minus Astro-curved family | 0.572 nat | `results/prl_reference_bridge/family_evidence.csv` |
| Hybrid3 best model gap: Cosmic superstrings minus best Astro-curved model | 0.663 nat | `results/prl_reference_bridge/hybrid3_bridge_lnZ.csv` |

The gaps above are the scale at which source-identification language must be
defended.

## Calibration and Systematic Terms

| effect | size | source / latest local check |
|---|---:|---|
| Official-density reproduction: maximum absolute residual across the five-configuration sweep | 0.694 nat | `theory/PRL_submission_package_manifest.md` |
| SIGW-delta prior-boundary extension to `log10_f_peak <= -3` | +0.426 nat | `theory/PRL_submission_package_manifest.md` |
| QMC/TI direct-vs-TI residual, maximum across checked rows | 0.0093 nat | `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json` |
| QMC/TI residual vs nested robust mean, SIGW-Gaussian | 0.118 nat | `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json` |
| QMC/TI residual vs nested robust mean, SIGW-delta | 0.063 nat | `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json` |
| QMC/TI residual vs nested robust mean, cosmic superstrings | 0.041 nat | `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json` |
| QMC/TI residual vs nested robust mean, environmental SMBHB fixed gamma | 0.066 nat | `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json` |
| Official-density extended rerun gap vs published SIGW-delta reference under PTArcade BHB-prior baseline | -0.137 nat | `results/T2_NG15yr/bayes_factors/ptarcade_ceffyl_density_extended.json` |
| Official-density extended rerun gap vs published superstring reference under PTArcade BHB-prior baseline | -0.491 nat | `results/T2_NG15yr/bayes_factors/ptarcade_ceffyl_density_extended.json` |

## Interpretation

The leading NG15 model gap and the hybrid3 family gap are both comparable to
the identified calibration and prior-boundary systematics.  They are much
larger than the independent QMC/TI direct-vs-TI integration residuals, so the
issue is not basic quadrature failure.  The safe interpretation is that current
public PTA products identify low-frequency spectral curvature more robustly
than a unique source class.

## Manuscript Use

The defensible sentence is:

> The remaining order-unity evidence gaps are comparable to the identified
> calibration and prior-boundary systematics, so they should not be interpreted
> as robust source discrimination.

Do not state that the data prefer a unique astrophysical or cosmological
origin on the basis of these gaps.
