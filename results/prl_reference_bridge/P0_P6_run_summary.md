# P0--P6 run summary

Generated: 2026-04-22 18:49:19

## Status

P0--P6 completed under `profile=production` with no failed rows in the
bridge manifest, main evidence rankings, stability table, or QMC/TI table.

The bridge route is technically viable: all 11 bridge models passed the
local PPTA, local NG15-HD, local EPTA, and official NG15-HD adapter smoke
tests.

## Main evidence result

### local3

The top local three-PTA posterior-summary ranking is:

| Rank | Model | Family | ln Z | Delta to best |
|---:|---|---|---:|---:|
| 1 | Cosmic-Superstrings | Cosmo-other | -62.594 +/- 0.109 | 0.000 |
| 2 | SIGW-Delta | SIGW-like | -62.711 +/- 0.110 | -0.117 |
| 3 | PBH-SIGW-Analytic | SIGW-like | -62.774 +/- 0.114 | -0.179 |
| 4 | SIGW-Gaussian | SIGW-like | -63.046 +/- 0.114 | -0.452 |
| 5 | SMBHB-Eccentric | Astro-curved | -63.215 +/- 0.121 | -0.620 |
| 6 | SMBHB-BrokenPL | Astro-curved | -63.344 +/- 0.121 | -0.750 |
| 7 | SMBHB-Env | Astro-curved | -63.527 +/- 0.123 | -0.933 |

### hybrid3

Replacing the local NANOGrav factor likelihood with the official NG15-HD
`ceffyl` density preserves the same qualitative top tier:

| Rank | Model | Family | ln Z | Delta to best |
|---:|---|---|---:|---:|
| 1 | Cosmic-Superstrings | Cosmo-other | -154.377 +/- 0.109 | 0.000 |
| 2 | SIGW-Delta | SIGW-like | -154.469 +/- 0.109 | -0.092 |
| 3 | PBH-SIGW-Analytic | SIGW-like | -154.550 +/- 0.114 | -0.173 |
| 4 | SIGW-Gaussian | SIGW-like | -154.850 +/- 0.115 | -0.472 |
| 5 | SMBHB-Eccentric | Astro-curved | -155.040 +/- 0.121 | -0.663 |
| 6 | SMBHB-BrokenPL | Astro-curved | -155.209 +/- 0.121 | -0.832 |
| 7 | SMBHB-Env | Astro-curved | -155.316 +/- 0.123 | -0.939 |

Immediate interpretation: curved-SMBHB controls enter the same
order-unity evidence tier as the leading cosmological templates. The result
supports a "shape identification, not robust source identification" PRL
line. It does not support a PBH discovery line.

## Calibration anchor

The NG15 local-to-official shift is approximately a common normalization
offset for this bridge implementation, with model-to-model variation much
smaller than the absolute offset:

- SMBHB-PowerLaw: `delta_cal = -91.808 +/- 0.171`
- PBH-SIGW-Analytic: `delta_cal = -91.787 +/- 0.153`
- SIGW-Delta: `delta_cal = -91.808 +/- 0.155`
- Cosmic-Superstrings: `delta_cal = -91.798 +/- 0.156`
- PhaseTransition-Total: `delta_cal = -91.856 +/- 0.175`

This validates the need for explicit scale separation, but in this bridge
set it does not create a large model-dependent reweighting.

## Family evidence

Hybrid3 family evidence with equal family mass and equal within-family
weights:

| Family | ln Z_F |
|---|---:|
| SIGW-like | -154.610 |
| Astro-curved | -155.182 |
| Cosmo-other | -155.279 |
| Astro-simple | -155.650 |

The SIGW-like family is highest, but the Astro-curved family is only
0.572 nat lower. Removing the leading curved row still gives
`Astro-curved = -155.261`, so the curved-family competitiveness is not
carried by a single row only.

## Gamma projection

The curved-SMBHB posterior samples project to effective slopes:

- `SMBHB-Eccentric`: gamma_eff medians 3.729 (PPTA), 3.807 (NG15), 4.022 (EPTA)
- `SMBHB-BrokenPL`: gamma_eff medians 3.781 (PPTA), 3.866 (NG15), 4.103 (EPTA)
- `SMBHB-Env`: gamma_eff medians 3.720 (PPTA), 3.814 (NG15), 4.087 (EPTA)

These projections move in the right direction for frequency-window
curvature, but they do not by themselves prove that the full PRD gamma
tension is astrophysical.

## Frequency-cut scan

In the diagnostic frequency-cut runs, PBH-SIGW-Analytic leads at
`f <= 5, 8, 12, 20 nHz`, while the all-bin diagnostic returns
Cosmic-Superstrings and SIGW-Delta essentially tied at the top. Curved SMBHB
models remain within about 0.4--1.2 nat depending on the cutoff.

This supports the claim that low-frequency bins drive much of the
discrimination, but the exact leader changes with bin inclusion.

## Robustness

Seed/nlive stability for the hybrid3 top tier:

- Cosmic-Superstrings: `-154.412` to `-154.348`
- SIGW-Delta: `-154.601` to `-154.428`
- PBH-SIGW-Analytic: `-154.571` to `-154.528`
- SIGW-Gaussian: `-154.982` to `-154.850`

QMC/TI cross-checks for the top promoted models:

- Cosmic-Superstrings: QMC/TI `lnZ = -154.429`
- SIGW-Delta: QMC/TI `lnZ = -154.604`
- PBH-SIGW-Analytic: QMC/TI `lnZ = -154.666`

The top-tier ranking is numerically stable at the order needed for the
merged-PRL planning discussion. The QMC/TI values are slightly lower than
the dynesty means but preserve the same order-unity non-discrimination
message.

## Outputs

- `results/prl_reference_bridge/bridge_model_manifest.csv`: exists
- `results/prl_reference_bridge/ng15_calibration_anchor.csv`: exists
- `results/prl_reference_bridge/local3_bridge_lnZ.csv`: exists
- `results/prl_reference_bridge/hybrid3_bridge_lnZ.csv`: exists
- `results/prl_reference_bridge/family_evidence.csv`: exists
- `results/prl_reference_bridge/gamma_projection.csv`: exists
- `results/prl_reference_bridge/frequency_cut_evidence.csv`: exists
- `results/prl_reference_bridge/robustness_budget.csv`: exists
- `results/prl_reference_bridge/ti_crosscheck_top_models.csv`: exists

## Next discussion gate

The next writing decision is whether to pivot the PRL title/abstract toward:

> Current public PTA data identify a low-frequency-curved broad spectral
> shape, but not a unique source class.

The numbers now support this stronger merged route, provided the manuscript
clearly labels `hybrid3` as a posterior-summary/official-anchor bridge, not
a full IPTA timing-residual likelihood.
