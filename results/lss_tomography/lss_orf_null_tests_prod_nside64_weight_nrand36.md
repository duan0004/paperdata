# DESI LSS ORF Template Null Tests

Generated: 2026-04-25T09:52:37

## Scope

- Input maps: DESI BGS_BRIGHT-21.5 selection-corrected overdensity maps.
- Projection: NANOGrav 15-year pulsar-pair antenna response geometry.
- Claim boundary: geometry/null-control products only; no timing-residual likelihood or source-identification evidence.

## Configuration

- Map tag: `prod_nside64_weight_nrand36`.
- HEALPix nside: `64`.
- Random-map nulls per bin: `128`.
- RNG seed: `20260425`.

## Template Summary

| bin | valid pixels | mask fraction | pair std | corr with HD | ORF vector |
|---|---:|---:|---:|---:|---|
| `all` | 13866 | 0.282104 | 0.000552231 | -0.297171 | `data/LSS/orf_templates/desi_bgs_bright21p5_all_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |
| `z_0p1_0p2` | 13850 | 0.281779 | 0.000553975 | 0.233365 | `data/LSS/orf_templates/desi_bgs_bright21p5_z_0p1_0p2_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |
| `z_0p2_0p3` | 13860 | 0.281982 | 0.000498425 | -0.277504 | `data/LSS/orf_templates/desi_bgs_bright21p5_z_0p2_0p3_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |
| `z_0p3_0p4` | 13861 | 0.282003 | 0.000583307 | -0.298310 | `data/LSS/orf_templates/desi_bgs_bright21p5_z_0p3_0p4_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |
| `z_0p1_0p4` | 13866 | 0.282104 | 0.000552231 | -0.297171 | `data/LSS/orf_templates/desi_bgs_bright21p5_z_0p1_0p4_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |

## Random-Map Nulls

| bin | real corr(HD) | null mean | null std | z | abs-percentile | real l2 z |
|---|---:|---:|---:|---:|---:|---:|
| `all` | -0.297171 | 0.016284 | 0.225901 | -1.388 | 0.805 | 0.918 |
| `z_0p1_0p2` | 0.233365 | -0.064877 | 0.255651 | 1.167 | 0.570 | 0.944 |
| `z_0p2_0p3` | -0.277504 | 0.016967 | 0.251767 | -1.170 | 0.664 | 0.503 |
| `z_0p3_0p4` | -0.298310 | -0.008838 | 0.240921 | -1.202 | 0.719 | 1.241 |
| `z_0p1_0p4` | -0.297171 | -0.014756 | 0.233647 | -1.209 | 0.750 | 0.948 |

## Inter-Template Correlations

| template_a | template_b | pair-vector Pearson r |
|---|---|---:|
| `all` | `z_0p1_0p2` | -0.019491 |
| `all` | `z_0p2_0p3` | 0.612738 |
| `all` | `z_0p3_0p4` | 0.908466 |
| `all` | `z_0p1_0p4` | 1.000000 |
| `z_0p1_0p2` | `z_0p2_0p3` | -0.369789 |
| `z_0p1_0p2` | `z_0p3_0p4` | -0.225217 |
| `z_0p1_0p2` | `z_0p1_0p4` | -0.019491 |
| `z_0p2_0p3` | `z_0p3_0p4` | 0.384474 |
| `z_0p2_0p3` | `z_0p1_0p4` | 0.612738 |
| `z_0p3_0p4` | `z_0p1_0p4` | 0.908466 |

## Interpretation Gate

- Passing this test means the DESI maps can be converted into stable pulsar-pair ORF vectors.
- It does not mean that PTA data prefer an LSS-correlated component.
- The next required test is an anisotropic ENTERPRISE likelihood with isotropic and random-map null injections.
