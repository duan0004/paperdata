# PRL Submission Human Gate

**Date**: 2026-04-22  
**Target**: Physical Review Letters  
**Current package**: `theory/paper_prl_submission.tex` +
`theory/prl_supplement.tex`

## Completed

- Author metadata inserted:
  - Ran DUAN
  - chenli science
  - rduan@chenli.science
- Main text and supplement compile with Tectonic:
  - `theory/pdf/revtex/paper_prl_submission.pdf`
  - `theory/pdf/revtex/prl_supplement.pdf`
- Reviewer-risk cleanup applied:
  - title/abstract now use curved-SMBHB controls rather than overclaiming a
    fully astrophysical model family;
  - submission-facing sources no longer use internal milestone labels or
    preregistration language;
  - Fig. 1 Panel A states the row-specific baseline convention;
  - supplement includes prior ranges, selection rule, and QMC/TI residual
    convention notes.
- Evidence cross-check added:
  - `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json`
  - `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.md`
- Acknowledgments section has been removed from the main draft rather than
  leaving a placeholder.  Add one back only if there is real grant,
  collaboration, or computing-resource text.

## Still Required From Author

1. **Public code-release URL**  
   Current source still contains `https://github.com/duan0004/paperdata`.  This must be
   replaced by a real public repository, archive, or DOI before submission.

2. **Acknowledgments decision**  
   Either keep acknowledgments omitted, or provide exact text for people,
   funding, grants, and computing resources.  Do not invent this text.

3. **Final upload choice**  
   Decide whether the APS package should use inline `thebibliography` as now,
   or whether to convert references to a `.bib` file.

## Submission Gate

The scientific package is PRL-plausible and the local compile gate now passes.
The package is still not submission-ready until the public code-release URL is
real.
