# PRL Systematic Evidence Envelope

**Generated**: 2026-04-22 09:53:41
**JSON**: `results/T2_NG15yr/bayes_factors/prl_H10_systematic_envelope.json`

## Cosmological Official-Density Reproduction

| Model | mean lnB | sample std | mean residual | max abs residual |
|---|---:|---:|---:|---:|
| SIGW-Gaussian | `+4.520` | `0.185` | `+0.477` | `0.694` |
| SIGW-delta | `+3.930` | `0.199` | `+0.146` | `0.387` |
| Cosmic superstrings | `+3.558` | `0.145` | `-0.271` | `0.491` |

## Curved SMBHB Control Family

| Model | class | mean lnB | sample std | QMC/TI residual |
|---|---|---:|---:|---:|
| `ecc_supp_fixed_gamma` | eccentricity_inspired_suppression | `+3.839` | `0.043` | `+0.056` |
| `broken_pl_fixed_gamma` | broken_powerlaw_curvature | `+3.656` | `0.154` | `+0.063` |
| `env_fixed_gamma` | environmental_turnover | `+3.453` | `0.062` | `+0.017` |
| `env_free_gamma` | environmental_turnover | `+3.360` | `0.202` |  |
| `env_broadbend_fixed_gamma` | environmental_turnover | `+2.992` | `0.142` |  |
| `env_lowbend_fixed_gamma` | environmental_turnover | `+2.666` | `0.103` |  |

## Budget Summary

| Component | Envelope / diagnostic |
|---|---|
| Official-density reproduction | mean residual max `0.477` nat; single-config max `0.694` nat |
| QMC/TI cross-check | main max residual `0.117` nat; curved-family max residual `0.063` nat |
| SIGW-delta prior boundary | `Delta lnB <= 0.426` nat over scanned upper bounds |
| Curved SMBHB family | tested-representative range `+3.453` to `+3.839`; family pass `True` |
| Frequency-bin driver | max `|Delta lnB(N=8 to 14)| = 0.687` nat; most models within 0.5 nat by `N=8` |

## Source-Identification Consequence

The robust cosmological means span `+3.558` to `+4.520`.  The best row in each curved-SMBHB class spans `+3.453` to `+3.839`.  The ranges overlap: `True`.

The top cosmological mean exceeds the top curved-SMBHB mean by `0.681` nat, while SIGW-delta exceeds the top curved-SMBHB row by only `0.091` nat. Therefore the PRL claim should be framed as non-discrimination on the current calibrated evidence scale, not as a new-physics source detection.
