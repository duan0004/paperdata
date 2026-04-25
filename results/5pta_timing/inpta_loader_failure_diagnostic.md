# InPTA Loader Failure Diagnostics

Generated: 2026-04-25T10:03:55

## Scope

- Read-only inspection of staged public InPTA DR1 par/tim files.
- No `.par` or `.tim` file is edited.
- Diagnostics identify likely loader issues; they do not define a production timing policy.

## Summary

- InPTA par/tim targets inspected: `28`.
- DMXR/DMXF fit-flag anomalies found: `1`.
- NB fit-flag anomalies found: `1`.
- WB tim files without `FORMAT 1`: `14`.
- Malformed `C`-prefix tim rows: `5`.

## Focused Load Attempts

| pulsar | band | FROM_PAR | DE440 | DE436 |
|---|---|---|---|---|
| `J0751+1807` | `NB` | `FAIL` | `FAIL` | `FAIL` |
| `J0751+1807` | `WB` | `FAIL` | `FAIL` | `FAIL` |
| `J0437-4715` | `NB` | `PASS` | `PASS` | `PASS` |
| `J0437-4715` | `WB` | `FAIL` | `FAIL` | `FAIL` |

## Interpretation

- `J0751+1807.NB.par` contains a fitted `DMXR1_0007` line; libstempo reports `param_dmxr1` during the failing full loader run, so this is the leading read-only explanation for the NB design-matrix failure.
- The WB `.tim` files lack the `FORMAT 1` header used by the NB files and all focused WB ENTERPRISE load attempts fail before usable TOAs are exposed.
- A production fix should either use an official InPTA loading/preprocessing recipe or define a documented exclusion/preprocessing policy; this diagnostic intentionally does not patch the timing files.

## Output Files

- `results/5pta_timing/inpta_par_fitflag_diagnostics.csv`
- `results/5pta_timing/inpta_tim_format_diagnostics.csv`
- `results/5pta_timing/inpta_tim_malformed_c_prefix.csv`
- `results/5pta_timing/inpta_focused_enterprise_loads.csv`
