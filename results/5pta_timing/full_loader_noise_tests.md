# Full Public-PTA Loader And Noise Tests

Generated: 2026-04-25T10:03:16

## Scope

- Read-only ENTERPRISE/enterprise_extensions loader gate.
- No timing model `.par` files are edited.
- This is not a timing-level evidence calculation.

## Loader Results

| test | status | expected pairs | loaded | failures | JSON |
|---|---|---:|---:|---:|---|
| `NG15` | `PASS` | 67 | 67 | 0 | `results/5pta_timing/full_loader_noise_NG15.json` |
| `EPTA_DR2new+` | `PASS` | 25 | 25 | 0 | `results/5pta_timing/full_loader_noise_EPTA_DR2new+.json` |
| `PPTA_DR3` | `PASS` | 32 | 32 | 0 | `results/5pta_timing/full_loader_noise_PPTA_DR3.json` |
| `InPTA_DR1_NB` | `PARTIAL` | 14 | 13 | 1 | `results/5pta_timing/full_loader_noise_InPTA_DR1_NB.json` |
| `InPTA_DR1_WB` | `FAIL` | 14 | 0 | 14 | `results/5pta_timing/full_loader_noise_InPTA_DR1_WB.json` |
| `MPTA_4p5yr` | `PASS` | 83 | 83 | 0 | `results/5pta_timing/full_loader_noise_MPTA_4p5yr.json` |

## Noise Product Inventory

| array | noise root | noise JSON files | auxiliary entries | note |
|---|---|---:|---:|---|
| `NG15` | `data/NG15yr/tutorials/data/15yr_wn_dict.json` | 1 |  | Official NANOGrav 15yr white-noise dictionary used by local NG15 loaders. |
| `EPTA_DR2new+` | `data/EPTA_DR2/epta-dr2/EPTA-DR2/noisefiles/DR2new+` | 25 | 3 | EPTA DR2new+ per-pulsar noise JSON plus red/dm/chrom dictionaries. |
| `PPTA_DR3` | `data/PPTA_DR3/ppta_dr3/toas_and_parameters/noisefiles` | 32 |  | PPTA DR3 single-pulsar noise JSON files. |
| `InPTA_DR1` | `data/InPTA_DR1/InPTA.DR1` | 0 | 28 | No explicit noise JSON was located in the staged public repository; DM time-series products are present. |
| `MPTA_4p5yr` | `data/MPTA_4p5yr/partim_datacentral_20260425/partim` | 0 |  | Staged public product is the noise-subtracted residual par/tim release; no separate noise JSON was located. |

## Gate Interpretation

- PASS means every discovered par/tim or feather product in that test loaded through the standard loader.
- PARTIAL means the public product is staged but at least one loader edge case must be resolved or excluded before production evidence.
- Arrays without explicit noise JSON require a documented likelihood policy before a full five-PTA evidence claim.
