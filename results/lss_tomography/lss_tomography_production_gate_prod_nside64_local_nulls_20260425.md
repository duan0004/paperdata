# LSS Tomography Production Gate

Generated: 2026-04-25T09:59:44

## Decision Boundary

- NG15+2MPZ published-null reproduction is the first science gate.
- DESI/WISExSCOS redshift-binned ORFs are production inputs only until that null is reproduced.
- No positive LSS-correlated signal should be trusted or reported from this gate.

Overall status: `BLOCKED_BEFORE_LSS_SCIENCE_CLAIM`.
2MPZ null reproduction status: `NOT_REPRODUCED_YET`.

## 2MPZ Download

| catalog | status | bytes | path | query/error |
|---|---|---:|---|---|
| 2MPZ | `exists` | 12001478 | `data/LSS/2MPZ/twompz_ra_dec_zphoto.csv.gz` | SELECT ra,dec,zPhoto FROM TWOMPZ..twompzPhotoz |

## Map Products

| catalog | bin | counts | valid fraction | delta std | delta map |
|---|---|---:|---:|---:|---|
| `twompz` | `all` | 933447 | 1.000000 | 0.611712 | `data/LSS/maps/twompz_all_prod_nside64_local_nulls_20260425_delta.npy` |
| `twompz` | `z_0_0p1` | 630986 | 1.000000 | 0.72704 | `data/LSS/maps/twompz_z_0_0p1_prod_nside64_local_nulls_20260425_delta.npy` |
| `twompz` | `z_0p1_0p2` | 279884 | 1.000000 | 0.751231 | `data/LSS/maps/twompz_z_0p1_0p2_prod_nside64_local_nulls_20260425_delta.npy` |
| `twompz` | `z_0_0p2` | 910870 | 1.000000 | 0.619961 | `data/LSS/maps/twompz_z_0_0p2_prod_nside64_local_nulls_20260425_delta.npy` |
| `twompz` | `z_0p2_0p3` | 21975 | 1.000000 | 1.6995 | `data/LSS/maps/twompz_z_0p2_0p3_prod_nside64_local_nulls_20260425_delta.npy` |
| `wisexscos` | `all` | 1.33761e+07 | 0.679830 | 0.658474 | `data/LSS/maps/wisexscos_all_prod_nside64_local_nulls_20260425_delta.npy` |
| `wisexscos` | `z_0_0p1` | 1.67762e+06 | 0.679830 | 0.85085 | `data/LSS/maps/wisexscos_z_0_0p1_prod_nside64_local_nulls_20260425_delta.npy` |
| `wisexscos` | `z_0p1_0p2` | 5.24673e+06 | 0.679830 | 0.689484 | `data/LSS/maps/wisexscos_z_0p1_0p2_prod_nside64_local_nulls_20260425_delta.npy` |
| `wisexscos` | `z_0_0p2` | 6.92436e+06 | 0.679830 | 0.695138 | `data/LSS/maps/wisexscos_z_0_0p2_prod_nside64_local_nulls_20260425_delta.npy` |
| `wisexscos` | `z_0p2_0p3` | 5.1625e+06 | 0.679830 | 0.680594 | `data/LSS/maps/wisexscos_z_0p2_0p3_prod_nside64_local_nulls_20260425_delta.npy` |

## NG15 ORF Templates

| catalog | bin | corr with HD | null z | abs-percentile | ORF vector |
|---|---|---:|---:|---:|---|
| `twompz` | `all` | 0.428995 | 2.109 | 0.977 | `data/LSS/orf_templates/twompz_all_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |
| `twompz` | `z_0_0p1` | 0.494065 | 2.281 | 1.000 | `data/LSS/orf_templates/twompz_z_0_0p1_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |
| `twompz` | `z_0p1_0p2` | 0.309055 | 1.393 | 0.719 | `data/LSS/orf_templates/twompz_z_0p1_0p2_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |
| `twompz` | `z_0_0p2` | 0.435749 | 2.089 | 0.953 | `data/LSS/orf_templates/twompz_z_0_0p2_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |
| `twompz` | `z_0p2_0p3` | 0.184503 | 0.891 | 0.523 | `data/LSS/orf_templates/twompz_z_0p2_0p3_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |
| `wisexscos` | `all` | 0.167460 | 0.618 | 0.406 | `data/LSS/orf_templates/wisexscos_all_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |
| `wisexscos` | `z_0_0p1` | -0.156699 | -0.628 | 0.383 | `data/LSS/orf_templates/wisexscos_z_0_0p1_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |
| `wisexscos` | `z_0p1_0p2` | 0.168852 | 0.781 | 0.445 | `data/LSS/orf_templates/wisexscos_z_0p1_0p2_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |
| `wisexscos` | `z_0_0p2` | 0.112917 | 0.584 | 0.344 | `data/LSS/orf_templates/wisexscos_z_0_0p2_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |
| `wisexscos` | `z_0p2_0p3` | 0.223224 | 0.967 | 0.594 | `data/LSS/orf_templates/wisexscos_z_0p2_0p3_prod_nside64_local_nulls_20260425_ng15_orf_pairs.npy` |

## DESI ORF Inventory

| status | pair count | pair std | ORF vector |
|---|---:|---:|---|
| `present` | 2211 | 0.000552231 | `data/LSS/orf_templates/desi_bgs_bright21p5_all_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |
| `present` | 2211 | 0.000553975 | `data/LSS/orf_templates/desi_bgs_bright21p5_z_0p1_0p2_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |
| `present` | 2211 | 0.000552231 | `data/LSS/orf_templates/desi_bgs_bright21p5_z_0p1_0p4_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |
| `present` | 2211 | 0.000498425 | `data/LSS/orf_templates/desi_bgs_bright21p5_z_0p2_0p3_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |
| `present` | 2211 | 0.000583307 | `data/LSS/orf_templates/desi_bgs_bright21p5_z_0p3_0p4_prod_nside64_weight_nrand36_ng15_orf_pairs.npy` |

## Blockers Before Science Use

1. Reproduce the published NG15+2MPZ null with the same 2MPZ template, selection/mask treatment, and likelihood definition.
2. Add isotropic and random-map null injections inside the ENTERPRISE anisotropy likelihood.
3. Run posterior sampling/evidence only after the null and injection gates pass.
