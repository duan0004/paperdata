# PRL Submission Human Gate

**Date**: 2026-04-23
**Target**: Physical Review Letters
**Current package**: `theory/paper_prl_submission.tex` +
`theory/prl_supplement.tex`

## Completed

- Author metadata inserted:
  - Ran DUAN
  - National Astronomical Observatories, Chinese Academy of Sciences
  - duanran@nao.cas.cn
- Main text and supplement compile with Tectonic:
  - `theory/pdf/revtex/paper_prl_submission.pdf`
  - `theory/pdf/revtex/prl_supplement.pdf`
- Reviewer-risk cleanup applied:
  - title/abstract now use curved-SMBHB controls rather than overclaiming a
    fully astrophysical model family;
  - submission-facing sources no longer use internal milestone labels or
    preregistration language;
  - Fig. 1 Panel A now uses only the matched PTArcade BHB-prior baseline;
  - supplement includes prior ranges, selection rule, and QMC/TI residual
    convention notes.
- Public GitHub code/data release inserted:
  - Zenodo version DOI `10.5281/zenodo.19688471`
  - Zenodo version URL `https://doi.org/10.5281/zenodo.19688471`
  - Zenodo concept DOI `10.5281/zenodo.19688471`
  - Zenodo concept URL `https://doi.org/10.5281/zenodo.19688471`
  - `https://github.com/duan0004/paperdata`
  - release tag `v1.0.4`
  - release URL:
    `https://github.com/duan0004/paperdata/releases/tag/v1.0.4`
- Evidence cross-check added:
  - `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json`
  - `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.md`
- Cross-PTA bridge-framing cleanup added:
  - `hybrid3` is now explicitly described as an anchored source-ranking
    stress test, not a fully calibrated multi-PTA Bayes-factor scale;
  - sequential bridge ablation is saved in
    `results/prl_reference_bridge/sequential_bridge_ablation.md`;
  - PTArcade/ceffyl version and plateau diagnostics are saved in
    `results/T2_NG15yr/bayes_factors/prl_ceffyl_plateau_diagnostic.md`;
  - final stability diagnostics are saved in
    `results/T2_NG15yr/bayes_factors/prl_final_stability_diagnostics.md`;
  - the bridge revision is mirrored at GitHub tag
    `prl-calibration-bridge-2026-04-23`.
- Acknowledgments section has been removed from the main draft rather than
  leaving a placeholder.  Add one back only if there is real grant,
  collaboration, or computing-resource text.

## Still Required From Author

1. **Acknowledgments decision**  
   Either keep acknowledgments omitted, or provide exact text for people,
   funding, grants, and computing resources.  Do not invent this text.

2. **Final upload choice**  
   Decide whether the APS package should use inline `thebibliography` as now,
   or whether to convert references to a `.bib` file.

## Submission Gate

The scientific package is PRL-plausible and the local compile gate now passes.
The Zenodo DOI, public GitHub code/data URL, and release URL are real.
Remaining human-level decision: whether acknowledgments are needed before
final APS upload.
