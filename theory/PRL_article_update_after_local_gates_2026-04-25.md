# PRL Article Update After Local Gates

Generated: 2026-04-25

## Scope

This update incorporates the locally runnable follow-up gates into the PRL
manuscript and supplement without upgrading any blocked diagnostic into a
science claim.

## Experiments Re-run Or Refreshed

- `python3 code/prl_final_stability_diagnostics.py`
  - Output: `results/T2_NG15yr/bayes_factors/prl_final_stability_diagnostics.md`
  - Confirms the curvature-mode compression, posterior-summary PPC, and
    evidence/systematics budget used in the current PRL framing.
- `python3 code/smbhb_population_controls.py`
  - Output: `theory/SMBHB_population_controls.md`
  - Reconfirms that the environmental turnover row is the physical SMBHB
    control, while broken-PL and eccentricity-inspired rows are
    phenomenological curvature controls.
- `python3 code/prl_package_static_gate.py`
  - Output: `results/T2_NG15yr/prl_package_static_gate.md`
  - Result: `structural_pass=True`, `submission_ready=True`.
- `python3 code/lss_tomography_manifest.py`
  - Result: manifest, mask, map, and ORF input checks pass locally.
- `python3 code/lss_2mpz_namaster_map_gate.py --nside 64 --tag prod_nside64_namaster_update_20260425 --lmax 12`
  - Result: blocked with `ModuleNotFoundError: No module named 'pymaster'`.
  - Interpretation: local environment blocker for exact NaMaster/MAP
    reconstruction; not a scientific result.
- `python3 code/lss_2mpz_lowell_reference_gate.py --nside 64 --tag prod_nside64_lowell_update_20260425 --nulls 128 --seed 20260425`
  - Output: `results/lss_tomography/lss_2mpz_lowell_reference_gate_prod_nside64_lowell_update_20260425.md`
  - Result: `LOWELL_FALLBACK_GATE_PASS`; this remains a deterministic fallback
    geometry diagnostic, not reproduction of the NG15+2MPZ published null.

## Manuscript Changes

- Updated `theory/paper_prl_submission.tex` in the cross-PTA bridge paragraph:
  - It now states that public NG15/EPTA/PPTA/MPTA timing inputs pass local
    loader checks, while InPTA, direct-combination, and noise-policy
    harmonization remain open.
  - It also states that the NG15+2MPZ LSS published-null baseline is not yet
    reproduced.
  - These diagnostics remain excluded from evidence claims.

## Supplement Changes

- Renamed S17 to `Timing-Level Production and Development Gates`.
- Added a full public-PTA loader/noise gate table:
  - NG15 67/67 PASS.
  - EPTA DR2new+ 25/25 PASS.
  - PPTA DR3 32/32 PASS.
  - MPTA 4.5 yr 83/83 PASS.
  - InPTA DR1 NB 13/14 PARTIAL.
  - InPTA DR1 WB 0/14 FAIL.
- Added exact-baseline manifest interpretation:
  - MPTA local ingestion now matches the 83-pulsar target.
  - InPTA remains policy-blocked.
  - Direct-combination audit finds 222 occurrences, 121 canonical groups,
    56 duplicate groups, two alias groups, and two canonical-name failures.
- Added mechanical ENTERPRISE timing-pilot table:
  - NG15 timing route, 8 pulsars, GW4: PASS.
  - NG15+EPTA+PPTA, 12 pulsars, GW4: PASS.
  - NG15+EPTA+PPTA, 12 pulsars, GW8/RN5: PASS.
- Added new S18 `LSS Tomography Production Gate`:
  - Records staged 2MPZ/DESI/WISE inputs.
  - Records NaMaster/pymaster blocker.
  - Records low-ell fallback gate and its claim boundary.
- Renumbered minimal local verification commands to S19 and added the refreshed
  P1-P4 commands.

## Compile And Static Gates

- Main text compiled with Tectonic:
  - `theory/pdf/revtex/paper_prl_submission.pdf`
  - 4 pages, 130265 bytes.
- Supplement compiled with Tectonic:
  - `theory/pdf/revtex/prl_supplement.pdf`
  - 11 pages, 119519 bytes.
- Latest preview copies updated:
  - `theory/pdf/paper_prl_submission_latest_preview.pdf`
  - `theory/pdf/prl_supplement_latest_preview.pdf`
- Log scan found no LaTeX errors, undefined control sequences, undefined
  references/citations, missing figures, or overfull boxes.
- Static package gate after the manuscript update:
  - `structural_pass=True`
  - `submission_ready=True`

## Claim Boundary

- Current PRL evidence claim remains the calibrated NG15 official-density
  comparison plus posterior-summary hybrid bridge.
- LSS tomography remains blocked before science use until the published
  NG15+2MPZ null is reproduced with the exact map construction and likelihood.
- Five-PTA timing-level family evidence remains blocked until exact public
  baseline reproduction, direct-combination policy, InPTA preprocessing, and
  noise-policy harmonization are all resolved.
