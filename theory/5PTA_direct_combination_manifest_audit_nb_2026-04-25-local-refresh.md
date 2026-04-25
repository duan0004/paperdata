# Five-PTA Direct-Combination Manifest Audit

Generated: 2026-04-25T10:03:39

## Scope

This is a read-only coordinate-grouping audit for the public five-PTA
direct-combination gate.  It does not edit `.par` files, does not write
converted timing models, and does not build a likelihood.

## Summary

- InPTA product used for this audit: `NB`.
- Array occurrences: `222`.
- Canonical pulsar groups: `121`.
- Duplicate canonical groups: `56`.
- Coordinate-canonical fallback failures: `2`.

The canonical group count matches the 121-pulsar public five-PTA target
once B/J aliases are grouped by coordinates.  This does not resolve the
published baseline policy; it only fixes the name-level manifest count.

## Alias Groups

| canonical name | input names | arrays | occurrences |
|---|---|---|---:|
| `J1857+0943` | `B1855+09;J1857+0943` | `EPTA;InPTA_NB;NG15;PPTA` | 4 |
| `J1939+2134` | `B1937+21;J1939+2134` | `InPTA_NB;NG15;PPTA` | 3 |

## Coordinate Fallbacks

| array | input name | fallback | error |
|---|---|---|---|
| `PPTA` | `J2241-5236` | `J2241-5236` | `ValueError: Some FBn parameters are set but FB0 is not.` |
| `MPTA` | `J1825-0319` | `J1825-0319` | `ValueError: Companion mass M2 cannot be negative (-0.4480919945782576 solMass)` |

## Outputs

- Occurrence CSV: `results/5pta_timing/direct_combination_manifest_audit_nb_2026-04-25-local-refresh_occurrences.csv`.
- Group CSV: `results/5pta_timing/direct_combination_manifest_audit_nb_2026-04-25-local-refresh_groups.csv`.
- JSON: `results/5pta_timing/direct_combination_manifest_audit_nb_2026-04-25-local-refresh.json`.

## Gate Decision

Status: `NAME_AND_COORDINATE_MANIFEST_READY_POLICY_BLOCKED`.

The local public timing inputs now reproduce the expected 121 canonical
pulsar groups at the manifest level.  Timing-level family evidence remains
blocked until the published direct-combination, InPTA preprocessing, and
noise-policy baseline is reproduced.
