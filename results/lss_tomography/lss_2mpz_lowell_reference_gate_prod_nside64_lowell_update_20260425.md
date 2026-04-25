# 2MPZ Low-Ell Reference-Template Fallback Gate

Generated: 2026-04-25T10:28:33

## Scope

- Target workflow: NG15+2MPZ published-null reproduction with 2MPZ slices `0<z<=0.1` and `0.1<z<=0.2`.
- Local implementation: Galactic latitude cut, apodized support, weighted monopole/dipole subtraction, low-ell truncation, and NG15 ORF projection.
- Claim boundary: this is a fallback template gate, not a reproduction of the published null and not LSS evidence.

## Blocking Difference From Reference

- `pymaster`/NaMaster is not installed locally, so exact C1 apodization and MASTER/MAP low-ell reconstruction are not available.
- This gate uses HEALPix mask smoothing and weighted `map2alm` truncation as a deterministic fallback.

## Configuration

- Tag: `prod_nside64_lowell_update_20260425`.
- HEALPix nside: `64`.
- Galactic cut: `|b| >= 20.0 deg`.
- Apodization fallback FWHM: `2.0 deg`.
- Support threshold: `0.001`.
- Low-ell range: `2 <= ell <= 12`.
- Random-map nulls per slice: `128`.

## 2MPZ Input

| rows read | rows valid | z min | z max | input |
|---:|---:|---:|---:|---|
| 934175 | 933447 | 7.57765e-06 | 0.401001 | `data/LSS/2MPZ/twompz_ra_dec_zphoto.csv.gz` |

## Map Products

| slice | counts | support fraction | delta rms | low-ell rms before norm | low-ell map |
|---|---:|---:|---:|---:|---|
| `z_0_0p1` | 630991 | 0.705363 | 0.647445 | 0.176256 | `data/LSS/maps/twompz_z_0_0p1_prod_nside64_lowell_update_20260425_lowell.npy` |
| `z_0p1_0p2` | 279879 | 0.705363 | 0.616377 | 0.117014 | `data/LSS/maps/twompz_z_0p1_0p2_prod_nside64_lowell_update_20260425_lowell.npy` |

## NG15 Pair ORFs

| slice | corr with HD | pair std | ORF vector |
|---|---:|---:|---|
| `z_0_0p1` | -0.102724 | 0.010563 | `data/LSS/orf_templates/twompz_z_0_0p1_prod_nside64_lowell_update_20260425_ng15_orf_pairs.npy` |
| `z_0p1_0p2` | -0.212425 | 0.0137923 | `data/LSS/orf_templates/twompz_z_0p1_0p2_prod_nside64_lowell_update_20260425_ng15_orf_pairs.npy` |

## Shuffle Nulls

| slice | real corr(HD) | null mean | null std | z | abs-percentile |
|---|---:|---:|---:|---:|---:|
| `z_0_0p1` | -0.102724 | -0.539892 | 0.141677 | 3.086 | 0.008 |
| `z_0p1_0p2` | -0.212425 | -0.243049 | 0.203512 | 0.150 | 0.336 |

## Decision Boundary

- Status: `LOWELL_FALLBACK_GATE_PASS`.
- The next gate is to run the fixed-parameter NG15 likelihood scan with these ORFs.
- DESI/WISExSCOS remain blocked from science interpretation until the exact NG15+2MPZ published null is reproduced.
