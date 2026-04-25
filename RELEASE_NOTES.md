# PRL Submission Code/Data Release v1.0.4

This release accompanies the PRL submission:

**Calibrated PTA Evidence: Low-Frequency Curvature Without a Unique Nanohertz Source**  
Ran DUAN, National Astronomical Observatories, Chinese Academy of Sciences, duanran@nao.cas.cn

It contains:

- PRL analysis scripts and static-gate checks.
- Official public NANOGrav/PTArcade input data products used by the analysis.
- Compact derived evidence tables, covariance diagnostics, QMC/TI checks, and figures.
- Posterior-summary bridge products, family-compression diagnostics, and final
  PRL stability checks.
- Non-evidence timing-level and LSS production-gate diagnostics used to define
  follow-up blockers.
- Manuscript and supplement sources plus compiled PDFs.
- `REPRODUCIBILITY.md`, `DATA_MANIFEST.md`, `environment.yml`, `.zenodo.json`, and `CITATION.cff`.

Large local MCMC chain dumps are intentionally excluded; summary diagnostics are included.

Expected local static gate:

```text
structural_pass=True
submission_ready=True
```

The timing-level and LSS-tomography gate files in this release are diagnostic
only. They do not enter any Bayes-factor or source-identification claim in the
Letter.
