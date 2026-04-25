# PRL Evidence Cross-Check: Sobol-QMC Thermodynamic Integration

**Generated**: 2026-04-25 09:33:45
**JSON**: `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json`

## Method

- Likelihood: official PTArcade `ceffyl` density.
- Models: PTArcade BHB-prior baseline, SIGW-Gaussian, SIGW-delta, cosmic superstrings, and fixed-gamma environmental SMBHB.
- QMC: `16384` Sobol points per scramble, `5` independent scrambles.
- TI grid: `81` beta values from `0` to `1`, with logarithmic spacing above `1.0e-04`.
- Evidence estimators: direct `log E_prior[L]` and thermodynamic integration `int_0^1 <log L>_beta d beta`.

## Main Comparison

| Model | QMC direct lnB | QMC TI lnB | dynesty robust mean | TI residual |
|---|---:|---:|---:|---:|
| SIGW-Gaussian | `+4.638 +/- 0.078` | `+4.637 +/- 0.078` | `+4.520` | `+0.117` |
| SIGW-delta | `+3.993 +/- 0.065` | `+3.992 +/- 0.065` | `+3.930` | `+0.063` |
| Cosmic superstrings | `+3.599 +/- 0.053` | `+3.598 +/- 0.053` | `+3.558` | `+0.040` |
| Environmental SMBHB, fixed gamma | `+3.519 +/- 0.082` | `+3.517 +/- 0.082` | `+3.453` | `+0.064` |

## Gate Interpretation

This cross-check is independent of dynesty and uses a power-posterior
thermodynamic-integration identity evaluated directly in the prior cube.
It is suitable for these 2--3 dimensional official-density integrals.
It is not a full parallel-tempering PTA-chain replacement.

Diagnostic pass criterion used here: the QMC-TI mean lnB for each main
model must agree with the existing dynesty robust mean to within `0.25` nat,
and the direct-QMC/TI discrepancy must remain below `0.05` nat.

**Diagnostic pass**: `True`

## Remaining Internal-Gate Note

The project AGENTS rule says Bayes factors must be cross-validated by
thermodynamic integration.  This file supplies a thermodynamic-integration
cross-check for the final official-density likelihood.  If a referee or
internal arbiter specifically requires thermodynamic integration from
parallel-tempered PTA chains, that remains a separate, larger run.
