# NANOGrav 15yr Data Download Instructions

## What Is Already Present (from git clone)

Repository cloned from: https://github.com/nanograv/15yr_stochastic_analysis

Contents:
- `tutorials/` — Tutorial notebooks and reduced data subset
  - `tutorials/data/15yr_wn_dict.json` — Official 15yr white noise dictionary (noise budget)
  - `tutorials/data/15yr_emp_distr.json` — Empirical distributions
  - `tutorials/data/feathers/` — 69 pulsar timing data files in Feather format (all 68 MSPs)
  - `tutorials/data/par/B1855+09_PINT_20220301.nb.par` — Example .par file
  - `tutorials/data/tim/B1855+09_PINT_20220301.nb.tim` — Example .tim file
- `data_release/` — Figure reproduction notebooks and data for the GWB paper
- 15 Jupyter notebooks

**Note**: The tutorial feather files contain all 68 pulsars' timing data but in reduced form. They are sufficient for running the tutorial notebooks and learning the analysis pipeline.

## Full Data Release (638.7 MB archive)

### Zenodo Record
- **DOI**: 10.5281/zenodo.7967584 (v1) / record 16051178 (v2.1.0)
- **URL**: https://zenodo.org/record/7967584
- **File**: `NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz` (638.7 MB)
- **MD5**: 557d42dd8486a5f8272d90dec9b228a8

### Download Command
```bash
curl -L --progress-bar \
  "https://zenodo.org/api/records/16051178/files/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz/content" \
  -o "./data/NG15yr/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz"
```

### Alternative: zenodo_get
```bash
pip install zenodo_get
zenodo_get 7967584 -o ./data/NG15yr/
```

### Extract After Download
```bash
# First peek at structure:
tar -tzf NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz | head -30

# Extract:
tar -xzf NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz \
  -C ./data/NG15yr/
```

### Archive Contents
- `.par` files — pulsar timing parameters for 68 millisecond pulsars
- `.tim` files — times of arrival (TOAs) for narrowband and wideband data
- Profile templates (FITS and pickle format)
- MCMC noise chains
- Clock correction files
- Correlation matrices

## References

- Data paper: Agazie et al. 2023, ApJL 951 L9, arXiv:2306.16213
- Noise analysis: arXiv:2306.16214
- GWB detection: arXiv:2306.16213

## Status

- [x] Tutorial repo cloned (GitHub): contains 69 feather files + noise dict
- [x] Full archive download in progress (Zenodo, 638.7 MB)
- [ ] Archive extraction pending download completion
