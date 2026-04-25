# SMBHB Population-Control Status

Generated: 2026-04-25T10:27:22

## Scope

P5 separates physically interpretable SMBHB controls from phenomenological curvature controls. The current PRL evidence supports a calibrated curved-SMBHB control statement, not a complete SMBHB population-synthesis marginalization.

## NG15 Official-Density Scale

| Model | physical status | lnB vs PTArcade BHB | scope |
|---|---|---:|---|
| `environmental turnover, fixed gamma` | `physical control` | +3.453 +/- 0.062 | SMBHB environmental-hardening spectral control |
| `environmental turnover, free gamma` | `physical control with extra spectral freedom` | +3.360 +/- 0.202 | environmental-hardening control plus free effective gamma |
| `environmental turnover, low-bend prior` | `physical-prior sensitivity` | +2.666 +/- 0.103 | environmental-hardening prior-range sensitivity |
| `environmental turnover, broad-bend prior` | `physical-prior sensitivity` | +2.992 +/- 0.142 | environmental-hardening prior-range sensitivity |
| `broken-power-law curvature` | `phenomenological curvature control` | +3.656 +/- 0.154 | low-frequency suppression surrogate, not population synthesis |
| `eccentricity-inspired suppression` | `phenomenological eccentricity-inspired control` | +3.839 +/- 0.043 | effective low-frequency suppression, not population synthesis |

## Bridge/Timing-Side Status

| Model | family | ln Z | delta to best | status |
|---|---|---:|---:|---|
| `SMBHB-Eccentric` | `Astro-curved` | -155.040 +/- 0.121 | -0.663 | bridge stress-test row; not full timing-level Bayes factor |
| `SMBHB-BrokenPL` | `Astro-curved` | -155.209 +/- 0.121 | -0.832 | bridge stress-test row; not full timing-level Bayes factor |
| `SMBHB-Env` | `Astro-curved` | -155.316 +/- 0.123 | -0.939 | bridge stress-test row; not full timing-level Bayes factor |
| `SMBHB-Turnover` | `Astro-simple` | -155.425 +/- 0.124 | -1.048 | bridge stress-test row; not full timing-level Bayes factor |
| `SMBHB-PowerLaw` | `Astro-simple` | -155.941 +/- 0.124 | -1.564 | bridge stress-test row; not full timing-level Bayes factor |

## Decision

- Main-text safe: environmental turnover as a physical SMBHB control; curved-SMBHB family as a tested-representative control family.
- Supplement safe: broken-power-law and eccentricity-inspired controls as low-frequency curvature surrogates.
- Not yet safe: claiming full astrophysical population-synthesis marginalization. That requires external population priors or simulation-derived template weights and a rerun through the official-density/timing-level evidence gates.
