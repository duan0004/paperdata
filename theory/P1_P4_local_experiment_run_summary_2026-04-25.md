# P1-P4 Local Experiment Run Summary

Generated: 2026-04-25T10:13:44

Boundary: this document summarizes local runnable gates only. LSS fixed-grid/profile outputs are diagnostics, and 5PTA timing pilots are mechanical gates; neither is used as science evidence or timing-level family Bayes factors.

## P1 PRL Stability Reinforcement

- Ran posterior-summary identifiability/PPC: `python3 code/pta_identifiability_ppc.py`; output reported `cut_rows=15`, `ppc_bins=294`, `ppc_summary=21`.
- Re-ran official-density evidence extension: `python3 code/bayes_factors_ptarcade_ceffyl_density_extended.py`.
  - SIGW-delta PTArcade lnB vs BHB-prior baseline: `+3.647 +/- 0.169`; published gap `-0.137` nat.
  - Cosmic superstrings lnB vs BHB-prior baseline: `+3.338 +/- 0.168`; published gap `-0.491` nat.
- Ran independent Sobol-QMC/TI cross-check: `python3 code/prl_evidence_ti_qmc_crosscheck.py`.
  - Diagnostic pass: `True`; SIGW-Gaussian TI residual vs dynesty robust mean: `+0.117` nat.
- Built uncertainty budget: `theory/PRL_evidence_uncertainty_budget_2026-04-25.md`.

## P2 LSS Local Gate

- Ran single-model NG15+2MPZ resource matrix for HD, low slice, high slice, and two-slice models. Full 67-pulsar single-model resource gates pass locally.
- Ran pairwise HyperModel thresholds. Low-basis `gw=4,rn=2` reaches full 67; reference-like `gw=14,rn=30` is stable at 24 pulsars and blocks at 25 in the local gate.
- Ran short sampler mechanics: pair-low `gw=14,rn=30,maxpsr=24,nsteps=100` PASS; this is mechanics only, not convergence/evidence.
- Ran 2MPZ fixed-grid diagnostics across epsilon/rn settings. Results are boundary/profile diagnostics only and remain below the published-null reproduction standard.
- Ran local LSS tomography gate: status `BLOCKED_BEFORE_LSS_SCIENCE_CLAIM`, 2MPZ null reproduction `NOT_REPRODUCED_YET`.
- Ran ORF geometry nulls with 128 nulls per bin. Correlations are geometry diagnostics and are not LSS evidence.

## P3 5PTA Manifest And Mechanical Gates

- Ran full loader/noise refresh: NG15 67/67 PASS, EPTA 25/25 PASS, PPTA 32/32 PASS, MPTA 83/83 PASS; InPTA NB is 13/14 PARTIAL and InPTA WB is 0/14 FAIL in the current read-only loader policy.
- Ran manifest v2 exact baseline gate: `BLOCKED_BEFORE_TIMING_LEVEL_FAMILY_EVIDENCE`.
- Ran NB and WB direct-combination manifest audits. Both yield 222 occurrences, 121 canonical groups, 56 duplicate groups, 2 alias groups, and 2 canonical failures.
- Ran InPTA loader failure diagnostics: one NB `DMXR1` fit-flag anomaly, WB files without `FORMAT 1`, and malformed `C`-prefix rows were identified read-only.
- Ran mechanical timing pilots: NG15 timing 8 pulsars gw4 PASS; three-array 4 each gw4 PASS; three-array 4 each gw8+rn5 PASS.
- Created duplicate-policy worklist for all 56 duplicate/alias groups.

## P4 Manuscript/Archive Sync

- Regenerated decisive evidence figure PDF/PNG/data.
- Regenerated bridge evidence figure PDF/PNG.
- Created reproducibility index with command provenance and SHA-256 checksums.

## Decision

- Current PRL line remains calibrated NG15 official-density evidence plus posterior-summary hybrid bridge.
- LSS production remains blocked before science claim until NG15+2MPZ published null is reproduced.
- 5PTA timing-level family evidence remains blocked until exact public baseline reproduction, direct-combination policy, InPTA policy, and noise policy are harmonized.
