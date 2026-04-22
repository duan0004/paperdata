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
- Public GitHub code/data release inserted:
  - `https://github.com/duan0004/paperdata`
  - commit `4d3972e11a24ffee026a8e2214eb12dd19fbc2ea`
- Evidence cross-check added:
  - `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json`
  - `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.md`
- Acknowledgments section has been removed from the main draft rather than
  leaving a placeholder.  Add one back only if there is real grant,
  collaboration, or computing-resource text.

## Still Required From Author

1. **Zenodo DOI decision**  
   The GitHub release is public and fixed by commit hash.  If the final
   submission requires a Zenodo DOI, archive the repository on Zenodo and add
   the DOI to the Data and Code Availability sentence.  Do not invent this DOI.

2. **Acknowledgments decision**  
   Either keep acknowledgments omitted, or provide exact text for people,
   funding, grants, and computing resources.  Do not invent this text.

3. **Final upload choice**  
   Decide whether the APS package should use inline `thebibliography` as now,
   or whether to convert references to a `.bib` file.

## Submission Gate

The scientific package is PRL-plausible and the local compile gate now passes.
The public GitHub code/data URL is real.  Remaining human-level decision:
whether to add a Zenodo DOI and acknowledgments before final APS upload.
