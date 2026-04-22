# PRL Submission Package Manifest

**Date**: 2026-04-22  
**Status**: reviewer-risk cleanup applied; submission blocked by public code URL.  
**Current title**: Curved SMBHB Spectra Compete with New-Physics
Explanations in PTA Source Identification

## Main Package

| Component | File | Status |
|---|---|---|
| REVTeX working draft | `theory/paper_prl_submission.tex` | reviewer-risk cleanup applied |
| Compressed Markdown source | `theory/paper_prl_compressed_draft.md` | reviewer-risk cleanup applied |
| Supplement Markdown draft | `theory/prl_supplement_draft.md` | local Markdown source ready |
| Supplement TeX draft | `theory/prl_supplement.tex` | Tectonic compile pass |
| Main figure | `results/T2_NG15yr/figures/prl_decisive_evidence_figure.pdf` | exists; static path check passes |
| Compiled PRL PDF | `theory/pdf/revtex/paper_prl_submission.pdf` | Tectonic compile pass |
| Compiled supplement PDF | `theory/pdf/revtex/prl_supplement.pdf` | Tectonic compile pass |
| Reproducibility manifest | `REPRODUCIBILITY.md` | ready except public repo URL |
| Environment manifest | `environment.yml` | ready |
| Formal compile report | `theory/T3.12_PRL_formal_compile_report.md` | done |
| Human submission gate | `theory/PRL_submission_human_gate.md` | done; code URL still required |
| Cover letter draft | `theory/PRL_cover_letter_draft.md` | ready except code URL |
| PRL assessment | `theory/T3.11_PRL_submission_assessment.md` | historical assessment; superseded by current cleanup |
| Static TeX check | `theory/T3.11_PRL_tex_static_check.md` | PASS; superseded by Tectonic compile report |
| Package static gate | `results/T2_NG15yr/prl_package_static_gate.md` | structural PASS; submission not ready |

## Main Numerical Claims

All Bayes factors below are relative to the PTArcade BHB-prior baseline and
use the official PTArcade `ceffyl` density.

| Claim | Value | Source |
|---|---:|---|
| SIGW-Gaussian reproduction | mean `lnB=+4.520`, residual `+0.477 nat` | `prl_official_density_robustness.json` |
| SIGW-delta reproduction | mean `lnB=+3.930`, residual `+0.146 nat` | `prl_official_density_robustness.json` |
| Superstring reproduction | mean `lnB=+3.558`, residual `-0.271 nat` | `prl_official_density_robustness.json` |
| Curved SMBHB control family | tested-representative range `lnB=+3.453` to `+3.839`; three classes pass | `prl_H7_astrophysical_family.json` |
| Top cosmology minus top curved SMBHB control | `0.681 nat`; SIGW-delta minus top curved control `0.091 nat` | `prl_H10_systematic_envelope.json` |
| TI/QMC evidence cross-check | four main non-baseline rows agree with dynesty robust means within `0.12 nat` | `prl_evidence_ti_qmc_crosscheck.json` |
| SIGW-delta prior-boundary shift | `+0.426 nat` when `log10_f_peak<=-3` | `prl_H4_H5_H6_experiments.json` |
| Low-frequency bin driver | most rows within `0.5 nat` of 14-bin evidence by first 8 bins; SIGW-Gaussian by 10 bins | `prl_H9_bin_driver_analysis.json` |
| CAR covariance diagnostic | effective rank `13.56/14`, max offdiag `0.227` | `car_null_calibration.json` |
| HD convergence | `Rhat(log10_A)=1.047`, `Rhat(gamma)=1.047` | `T2_multichain_Rhat.json` |

## Supporting Material Candidates

The current supplement draft is `theory/prl_supplement_draft.md`.  It collects:

- local KDE failure and ablation table:
  `results/T2_NG15yr/bayes_factors/prl_ablation_table.md`;
- full 15-model official-density ranking:
  `results/T2_NG15yr/bayes_factors/prl_H4_H5_H6_experiments.md`;
- SIGW-delta prior-boundary scan:
  `results/T2_NG15yr/bayes_factors/prl_H4_H5_H6_experiments.md`;
- curved-SMBHB control family:
  `results/T2_NG15yr/bayes_factors/prl_H7_astrophysical_family.md`;
- official-density cumulative-bin driver:
  `results/T2_NG15yr/bayes_factors/prl_H9_bin_driver_analysis.md`;
- systematic evidence envelope:
  `results/T2_NG15yr/bayes_factors/prl_H10_systematic_envelope.md`;
- Sobol-QMC thermodynamic-integration evidence cross-check:
  `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.md`;
- CAR null calibration:
  `results/T2_NG15yr/covariance/car_null_calibration.md`;
- legacy local-KDE `lnB(N_bins)` truncation diagnostic:
  `results/T2_NG15yr/bayes_factors/Nbins_scan.json`;
- full template derivations:
  `theory/paper_full_draft.md` and source-specific template notes.

## Current Static Gate

Run:

```bash
python3 code/prl_package_static_gate.py
```

Latest output:

- `results/T2_NG15yr/prl_package_static_gate.json`
- `results/T2_NG15yr/prl_package_static_gate.md`

Latest verdict:

- structural pass: `True`;
- submission ready: `False`;
- blockers: public code-release URL unresolved;
- word counts: REVTeX source `1612`, compressed Markdown `1512`, supplement
  Markdown `986`, supplement TeX `1604`;
- background process scan: `0` matching compute tasks.

## Hard Blockers

These must not be invented:

1. public code-release URL;
2. acknowledgments, grant, or computing-resource text if the author wants them.
   The current PRL source omits acknowledgments rather than using a placeholder.

Technical blockers:

1. final PRL/arXiv package needs a PDF visual check after inserting the public
   code URL;
2. final arbiter should be run on the actual compiled package, not only the
   Markdown draft.

## Recommended Final Gate

Before submission:

1. replace `https://github.com/duan0004/paperdata`;
2. recompile `theory/paper_prl_submission.tex` and `theory/prl_supplement.tex`;
3. inspect the PDF table for two-column overflow;
4. confirm Figure 1 renders and is readable at PRL column width;
5. run the static citation/path check again;
6. run the final arbiter on the compiled package.
