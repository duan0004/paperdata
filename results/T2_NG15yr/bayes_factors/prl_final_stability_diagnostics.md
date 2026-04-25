# PRL Final Stability Diagnostics

Generated: 2026-04-25T10:27:22

These tables consolidate existing completed runs. They do not introduce new timing-level or nested-sampling evidence.

## Curvature-Mode Compression

| diagnostic | value | interpretation |
|---|---:|---|
| Astro-curved family gap on hybrid3 | `-0.572` nat | source-family separation remains order-unity |
| curved-SMBHB gamma_eff median range | `3.720`--`4.103` | window-dependent low-frequency curvature is directly visible |
| top model gap, NG15 official scale | `0.681` nat | comparable to calibration/prior systematics |

## Posterior-Summary PPC

| family | mean low-4 | max low-4 | mean high-bin | interpretation |
|---|---:|---:|---:|---|
| Astro-curved | `0.462` | `1.617` | `0.923` | posterior-summary diagnostic only; comparable low-frequency residuals across top families |
| Cosmo-other | `0.478` | `1.559` | `0.911` | posterior-summary diagnostic only; comparable low-frequency residuals across top families |
| SIGW-like | `0.465` | `1.646` | `0.914` | posterior-summary diagnostic only; comparable low-frequency residuals across top families |

## Evidence/Systematics Budget

| effect | size | scope |
|---|---:|---|
| official-density reproduction residual | `0.694` nat | maximum absolute single-configuration residual among the three calibrated NANOGrav targets |
| SIGW-delta peak-frequency prior boundary | `0.426` nat | scan over tested upper peak-frequency bounds |
| dynesty/QMC-TI residual, main rows | `0.117` nat | independent Sobol-QMC/TI cross-check |
| dynesty/QMC-TI residual, curved rows | `0.063` nat | family cross-check for tested curved-SMBHB representatives |
| hybrid3 leading-row seed/live range | `0.173` nat | maximum range over the bridge robustness sweep |
| SIGW-Gaussian minus best curved-SMBHB model | `0.681` nat | matched NG15 official-density scale |
| SIGW-like minus Astro-curved family | `0.572` nat | hybrid3 equal-weight tested-representative family mixture |

## Environmental Prior Sensitivity

| model | lnB | gap to SIGW-Gaussian | gap to SIGW-delta | gap to superstrings | status |
|---|---:|---:|---:|---:|---|
| `env_fixed_gamma` | `3.453` | `1.067` | `0.476` | `0.104` | physical control |
| `env_free_gamma` | `3.360` | `1.160` | `0.570` | `0.197` | physical control with extra spectral freedom |
| `env_broadbend_fixed_gamma` | `2.992` | `1.527` | `0.937` | `0.565` | physical-prior sensitivity |
| `env_lowbend_fixed_gamma` | `2.666` | `1.854` | `1.264` | `0.892` | physical-prior sensitivity |

## Interpretation

The final stability diagnostics support the PRL framing without upgrading the claim to full timing-level source identification. The public posterior-summary PPCs show comparable low-frequency residuals for SIGW-like, Astro-curved, and Cosmo-other leading families; the matched-scale evidence gaps remain comparable to identified calibration, prior-boundary, and bridge robustness systematics. Environmental-turnover prior variants remain competitive with at least part of the leading cosmological tier, but they do not constitute a complete population-synthesis marginalization.
