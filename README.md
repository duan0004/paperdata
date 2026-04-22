# Curved SMBHB Spectra Compete with New-Physics Explanations in PTA Source Identification

This repository contains the code, data products, derived evidence tables, and
manuscript files for the PRL submission:

**Curved SMBHB Spectra Compete with New-Physics Explanations in PTA Source
Identification**  
Ran DUAN, National Astronomical Observatories, Chinese Academy of Sciences, duanran@nao.cas.cn

## Contents

- `code/`: scripts used for the spectral-template checks, official-density
  evidence calculations, QMC/TI cross-checks, covariance diagnostics, and PRL
  figure generation.
- `data/NG15yr/`: public NANOGrav/PTArcade inputs used by the analysis,
  including the official PTArcade model files, official `ceffyl` density
  product, NANOGrav 15-year tutorial data, and pre-sampled cores.
- `results/T2_NG15yr/`: compact derived outputs used in the paper tables and
  figures.  Large local MCMC chain directories are intentionally not included;
  their summary JSON files are included.
- `theory/`: REVTeX manuscript source, supplement source, cover letter draft,
  and compiled PDFs.
- `environment.yml`: pinned Python environment used for the saved results.
- `REPRODUCIBILITY.md`: commands, input fingerprints, and expected outputs.
- `.zenodo.json` and `CITATION.cff`: metadata for the Zenodo/GitHub archived
  release.

## Public Data Provenance

The analysis uses:

- NANOGrav 15-year data set: arXiv:2306.16213.
- NANOGrav 15-year noise analysis: arXiv:2306.16214.
- NANOGrav new-physics comparison: arXiv:2306.16219.
- Official NANOGrav/PTArcade model files: Zenodo DOI
  `10.5281/zenodo.8084351`.
- Official PTArcade `ceffyl` density product: Zenodo DOI
  `10.5281/zenodo.10495907`.

Third-party data retain their original licenses and citation requirements.

## Minimal Verification

Create the environment:

```bash
conda env create -f environment.yml
conda activate pta-gwb-ng15
```

Run the fast checks:

```bash
python3 code/gwb_templates.py
python3 -m py_compile code/prl_*.py
python3 code/prl_decisive_evidence_figure.py
python3 code/prl_H10_systematic_envelope.py
python3 code/prl_package_static_gate.py
```

The expected static gate status is `structural_pass=True`.  The submission gate
is ready only when the manuscript points to this public repository.

## Main Reproducibility Targets

- Main figure:
  `results/T2_NG15yr/figures/prl_decisive_evidence_figure.pdf`
- Main evidence table:
  `results/T2_NG15yr/bayes_factors/prl_official_density_robustness.json`
  and `results/T2_NG15yr/bayes_factors/prl_H7_astrophysical_family.json`
- QMC/TI cross-check:
  `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json`
- Cumulative-bin scan:
  `results/T2_NG15yr/bayes_factors/prl_H9_bin_driver_analysis.json`
- Static package gate:
  `results/T2_NG15yr/prl_package_static_gate.json`

The full nested-sampling sweeps are reproducible from the included scripts and
data, but they are computationally heavier than the figure/static-gate checks.
