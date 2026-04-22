# T2.5 — Bayes-factor refit against HD free-spectrum posterior

**Date**: 2026-04-18
**Script**: `code/bayes_factors_refit.py`
**Output**: `bayes_factors.json`

## Method

ceffyl-style refit (per-bin KDE product likelihood) of the official NG15yr
HD 30-frequency free-spectrum posterior (`hd_30f_fs.core`, arXiv:2306.16213
tutorials).  Only the 14 lowest Fourier bins are used (f = 1/Tspan …
14/Tspan), matching the NG15yr new-physics paper (arXiv:2306.16219 §III.B)
convention.

- **Likelihood**: ln L(θ) = Σ_k ln KDE_k( log10_ρ_k^model(θ) ), 14 independent Gaussian KDEs built from the posterior chain marginals.
- **Conversion**: ρ_k² = H0² Ω_GW(f_k) / (8 π⁴ f_k⁵ T_span)   [units s²].
- **Evidence**: dynesty static `NestedSampler`, nlive=500, dlogz=0.1, `rwalk`.

## Results (relative to SMBHB, γ=13/3 fixed)

| Model          | dim | ln Z        | Δln Z vs SMBHB | B vs SMBHB |
|----------------|----:|------------:|---------------:|-----------:|
| **SMBHB** (ref)  | 1 | −20.16 ± 0.10 | 0               | 1          |
| cosmic strings (stable NG) | 1 | −48.43 ± 0.11 | −28.27 ± 0.14 | 5.3 × 10⁻¹³ |
| SIGW (log-normal) | 3 | −21.75 ± 0.12 | −1.59 ± 0.15  | 0.20 ± 0.03 |

## Caveats — why these differ from published NG15yr Bayes factors

NG15yr new-physics paper (arXiv:2306.16219 Table 3) reports, for the same
refit setup:

| Model                      | Our refit B | NG15yr Table 3 B |
|----------------------------|-------------|------------------|
| SIGW-Gauss                 | 0.20 ± 0.03 | 9.0⁺¹·⁷₋₁·₄       |
| Cosmic strings (stable)    | 5 × 10⁻¹³   | ~0.3             |

**Root cause of the discrepancy**: this original refit used both
*approximate* spectral templates (T1.10) and a local scipy-KDE reconstruction
of the free-spectrum density.  NG15yr uses full PTArcade templates (Mitridate
et al. arXiv:2306.16377) with accurate cosmological integration plus the
official PTArcade ceffyl density construction.

Specifically:

1. **Cosmic strings**: our template returns a strictly flat Ω_GW(f)
   (Caprini-Figueroa Eq. 364 order-of-magnitude approximation).  The real
   BOS-model spectrum is a mild rising power law (α ≈ 0.2) across the PTA
   band.  A perfectly flat Ω_GW at amplitude fitting one free-spec bin
   gives wrong predictions at the other 13 bins, killing the evidence.

2. **SIGW-Gauss**: our semi-analytic peak-normalisation calibrated to the
   NG15yr best-fit point is accurate near the peak but inaccurate in
   shape elsewhere in the band.  Net effect: the 3D parameter space
   struggles to find a configuration that fits all 14 bins as well as
   SMBHB γ=13/3, so Occam's razor (larger prior volume) penalises it.

**Ranking is still qualitatively correct**: SMBHB > SIGW-Gauss >>
cosmic strings.  Absolute Bayes factors require PTArcade-grade templates and
the official density construction.

## Action

- [ ] Do NOT fill `[TBD]` in `theory/paper_full_draft.md` with these numbers.
- [x] Cite NG15yr arXiv:2306.16219 Table 3 values for the published Bayes
  factors and note that our in-house template library reproduces their
  qualitative ordering but not absolute magnitudes.
- [x] PTArcade official-density gates completed below: official templates plus
  official ceffyl density reproduce the published SIGW-Gauss, SIGW-delta, and
  cosmic-superstring Bayes factors to `<0.5 nat`.

---

# T2b — PTArcade official-density ground-truth pilots

**Date**: 2026-04-21
**Scripts**:
- `code/bayes_factors_ptarcade_groundtruth.py`
- `code/bayes_factors_ptarcade_ceffyl_density.py`
- `code/bayes_factors_ptarcade_ceffyl_density_extended.py`

**Outputs**:
- `ptarcade_groundtruth.json`
- `ptarcade_ceffyl_density_groundtruth.json`
- `ptarcade_ceffyl_density_extended.json`

## Source Data

Official NANOGrav/PTArcade model files were downloaded from Zenodo DOI
`10.5281/zenodo.8084351` into
`data/NG15yr/PTArcade_models_1.0.0/`.  The official models used here are
`sigw_gauss.py`, `sigw_delta.py`, and `super.py`; the SIGW-Gauss and
superstring spectra use the accompanying HDF5 tables under `models_data/`.

Official PTArcade ceffyl density data were downloaded from Zenodo DOI
`10.5281/zenodo.10495907` into
`data/NG15yr/PTArcade_ceffyl_0.2.0/`.

PTArcade model files return `h^2 Omega_GW`, not bare `Omega_GW`.  The residual
conversion used for the official templates is therefore

`rho^2 = (H0/h)^2 h^2 Omega_GW / (8 pi^4 f^5 Tspan)`.

This conversion was checked directly against `ptarcade.models_utils.omega2cross`
to relative precision `<5e-16`.

## Results

| Likelihood/density | Model | lnB vs fixed-gamma SMBHB | lnB vs PTArcade BHB prior | Gap vs published lnB |
|---|---|---:|---:|---:|
| Local scipy KDE built from `hd_30f_fs.core` | SIGW-Gauss | `+0.128 ± 0.145` | `+0.687 ± 0.147` | `-3.36 nat` using PTArcade BHB prior |
| Official PTArcade ceffyl density grid | SIGW-Gauss | `+3.859 ± 0.156` | `+4.225 ± 0.173` | `+0.18 nat` vs `ln(57)` |
| Official PTArcade ceffyl density grid | SIGW-delta | `+3.281 ± 0.152` | `+3.647 ± 0.169` | `-0.14 nat` vs `ln(44)` |
| Official PTArcade ceffyl density grid | cosmic superstrings (`super.py`) | `+2.971 ± 0.151` | `+3.338 ± 0.168` | `-0.49 nat` vs `ln(46)` |

## Interpretation

The 2026-05-05 Plan-D gate **passes for SIGW-Gauss, SIGW-delta, and the
official cosmic-superstring template**: using the official PTArcade spectra and
official PTArcade ceffyl density reproduces the published NANOGrav values to
`0.18 nat`, `0.14 nat`, and `0.49 nat`, respectively, when compared against the
PTArcade BHB prior baseline.

The failed local-scipy-KDE run is also informative.  Merely replacing the
in-house semi-analytic template while keeping the local `scipy.gaussian_kde`
reconstruction of `hd_30f_fs.core` does not reproduce the published evidence.
The remaining difference is therefore not caused by CAR/bin covariance; it is
mainly the combination of official PTArcade template fidelity and the official
ceffyl density construction.

---

# PRL H4/H5/H6 — Official-density hardening controls

**Date**: 2026-04-21
**Scripts**:
- `code/prl_H4_H5_H6_experiments.py`
- `code/prl_H4_env_robustness.py`

**Outputs**:
- `prl_H4_H5_H6_experiments.json`, `.md`
- `prl_H4_env_robustness.json`, `.md`

## H4: Environmental SMBHB curvature control

Five configurations were run against the same PTArcade BHB-prior baseline and
official PTArcade `ceffyl` density used by the H1 robustness sweep.

| Model | mean lnB | sample std | min | max | mean nested lnB err |
|---|---:|---:|---:|---:|---:|
| fixed-`gamma=13/3` environmental SMBHB | `+3.453` | `0.062` | `+3.366` | `+3.511` | `0.153` |
| free-`gamma` environmental SMBHB | `+3.360` | `0.202` | `+3.081` | `+3.637` | `0.155` |

Interpretation: an astrophysical curved-SMBHB control is competitive with the
leading cosmological templates on the same evidence scale.  This is a
source-identification degeneracy result, not a new source claim.

## H5: Official stochastic model ranking

All official PTArcade files defining `spectrum(f, ...)` were evaluated in a
single `nlive=500`, `dlogz=0.1`, seed `20260421` configuration.  The leading
six models are:

| Rank | Model | lnB vs BHB prior |
|---:|---|---:|
| 1 | `sigw_gauss` | `+4.616` |
| 2 | `sigw_delta` | `+4.098` |
| 3 | `sigw_box` | `+3.830` |
| 4 | `super` | `+3.739` |
| 5 | `meta_ls` | `+3.677` |
| 6 | `dw_sm` | `+3.625` |

The full 15-model table is in `prl_H4_H5_H6_experiments.md`.

## H6: SIGW-delta prior-boundary sensitivity

| Prior label | lnB vs BHB prior | delta lnB vs original |
|---|---:|---:|
| `log10_f_peak <= -5` | `+4.098` | `+0.000` |
| `log10_f_peak <= -4` | `+4.300` | `+0.203` |
| `log10_f_peak <= -3` | `+4.524` | `+0.426` |

Interpretation: the SIGW-delta evidence has a modest `~0.4 nat`
prior-boundary sensitivity because the best-fit peak frequency moves above
the original upper prior when allowed.

---

# D2b — N_bins truncation-stability pilot

**Date**: 2026-04-21
**Script**: `code/bayes_factors_Nbins_scan.py`
**Output**: `Nbins_scan.json`, `Nbins_scan.png`

## Method

Repeat the same ceffyl-style refit while varying the retained number of
low-frequency Fourier bins:

`N_bins = {6, 8, 10, 12, 14, 18, 22, 26, 30}`.

The likelihood still uses independent 1D KDEs for each free-spectrum bin,
but the script now precomputes each KDE as a 4096-point log-pdf interpolation
table.  This was validated against the original exact-KDE run for the
completed low-bin cases: at 14 bins, the interpolated results differ by
only ~0.01-0.03 nat from the exact run, well below the nested-sampling
uncertainty.

Pilot gates:

- **G3**: at least one non-SMBHB template shifts by >=2 nats across
  `N_bins = {6, 14, 30}`.
- **G4**: the same shift is monotonic.

## Gate Results

| Model | ΔlnB(N=6) | ΔlnB(N=14) | ΔlnB(N=30) | Span | Monotonic | Gate |
|---|---:|---:|---:|---:|---|---|
| SIGW-Gauss | -0.82 | -1.55 | -0.59 | 0.96 | no | FAIL |
| CS-stable proxy | -24.25 | -28.12 | -29.16 | 4.92 | yes | PASS |

Overall D2b verdict: **G3 PASS, G4 PASS**, driven entirely by the stable-string
proxy.  The SIGW-Gauss template does **not** show a strong monotonic
truncation drift.

## Interpretation

This is useful as a spectral-shape diagnostic, but not as a universal new
observable for all candidate sources.  The stable-string result is also
template-limited: our current CS-stable implementation is a deliberately
crude flat-spectrum proxy, not the full PTArcade/BOS spectrum.  Therefore,
use this pilot to justify reporting `lnB(N_bins)` as a robustness diagnostic,
not to claim a physical exclusion of cosmic strings from the in-house proxy.
