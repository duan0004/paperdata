# Reproducibility Manifest

**Project**: PTA stochastic-GWB spectral-template comparison with NANOGrav 15yr data and cross-PTA posterior-summary bridge  
**Last updated**: 2026-04-23  
**Archived release DOI**: `10.5281/zenodo.19688587`  
**Archived release URL**: `https://doi.org/10.5281/zenodo.19688587`  
**Public repository**: `https://github.com/duan0004/paperdata`  
**GitHub release tag**: `v1.0.1`  
**Submission package commit**: `2bc1dfeee7d452f6c3b58a8fd8f339e6ff0091f6` for the archived NANOGrav-calibration snapshot; the P0--P6 bridge update is recorded in the Git history after that release.  
**Runtime used for saved results**: Python 3.9.6 on macOS 26.0 arm64

This manifest records the local environment, input data fingerprints, and
minimal commands needed to reproduce the saved validation and evidence results.

## Environment

The pinned environment file is `environment.yml`.

Observed package versions in the local runtime:

| Package | Version |
|---|---:|
| `enterprise-pulsar` | `3.4.4` |
| `enterprise_extensions` | `3.0.3` |
| `PTMCMCSampler` | `2.1.4` |
| `dynesty` | `3.0.0` |
| `ceffyl` | `1.41.2` |
| `ptarcade` | `1.1.5` |
| `la_forge` | `1.1.0` |
| `numpy` | `1.26.4` |
| `scipy` | `1.13.1` |
| `h5py` | `3.14.0` |
| `matplotlib` | `3.9.4` |
| `astropy` | `6.0.1` |
| `pandas` | `2.3.1` |
| `emcee` | `3.1.6` |
| `mpmath` | `1.3.0` |
| `sympy` | `1.14.0` |

## Data Inputs

| Input | Path | SHA256 |
|---|---|---|
| NANOGrav 15yr white-noise dictionary | `data/NG15yr/tutorials/data/15yr_wn_dict.json` | `85da2a6e8c727b30a502e8247378273b4ffd782e32103716ba8791e7481828ca` |
| Official NANOGrav/PTArcade model archive | `data/NG15yr/PTArcade_models_1.0.0/models_1.0.0.tar.gz` | `5354c4f766804f93711091161a2b26209c6abc75a0867350b3607ff196ae3845` |
| Official PTArcade ceffyl density archive | `data/NG15yr/PTArcade_ceffyl_0.2.0/ng15_30f_fs{hd}_ceffyl.zip` | `1a4d7ffba73eecedd19ffccbb9645141f50d38fea0e2a32c51d344461034a438` |

Additional NANOGrav pre-sampled cores used by the analysis are under
`data/NG15yr/tutorials/presampled_cores/`, especially `hd_30f_fs.core`.

## Saved Results

| Result | File |
|---|---|
| Warm HD production-chain summary | `results/T2_NG15yr/T2_summary_prodHD_final.json` |
| Cold HD chain summary | `results/T2_NG15yr/T2_summary_prodHD_cold_final.json` |
| Warm+cold Gelman-Rubin diagnostics | `results/T2_NG15yr/T2_multichain_Rhat.json` |
| Main local ceffyl-style refit | `results/T2_NG15yr/bayes_factors/bayes_factors.json` |
| KDE bandwidth / sampler sensitivity sweep | `results/T2_NG15yr/bayes_factors/sensitivity.json` |
| SMBHB gamma-prior-volume check | `results/T2_NG15yr/bayes_factors/smbhb_gamma_prior_volume.json` |
| `lnB(N_bins)` truncation-stability diagnostic | `results/T2_NG15yr/bayes_factors/Nbins_scan.json` |
| PTArcade SIGW-Gauss official-density gate | `results/T2_NG15yr/bayes_factors/ptarcade_ceffyl_density_groundtruth.json` |
| PTArcade SIGW-delta / superstring official-density gate | `results/T2_NG15yr/bayes_factors/ptarcade_ceffyl_density_extended.json` |
| PRL official-density robustness sweep | `results/T2_NG15yr/bayes_factors/prl_official_density_robustness.json` |
| PRL H4 environmental-SMBHB robustness | `results/T2_NG15yr/bayes_factors/prl_H4_env_robustness.json` |
| PRL H5/H6 official model/prior scan | `results/T2_NG15yr/bayes_factors/prl_H4_H5_H6_experiments.json` |
| PRL H7 astrophysical curvature family | `results/T2_NG15yr/bayes_factors/prl_H7_astrophysical_family.json` |
| PRL H8 decisive evidence figure data | `results/T2_NG15yr/figures/prl_decisive_evidence_figure_data.json` |
| PRL H9 cumulative-bin driver | `results/T2_NG15yr/bayes_factors/prl_H9_bin_driver_analysis.json` |
| PRL H10 systematic evidence envelope | `results/T2_NG15yr/bayes_factors/prl_H10_systematic_envelope.json` |
| PRL Sobol-QMC thermodynamic-integration evidence cross-check | `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json` |
| CAR covariance pilot | `results/T2_NG15yr/covariance/pilot_report.json` |
| CAR template projection diagnostic | `results/T2_NG15yr/covariance/template_projection.json` |
| CAR null calibration | `results/T2_NG15yr/covariance/car_null_calibration.json` |
| PRL package static gate | `results/T2_NG15yr/prl_package_static_gate.json` |
| P0--P6 bridge run summary | `results/prl_reference_bridge/P0_P6_run_summary.md` |
| local3 versus hybrid3 bridge ranking | `results/prl_reference_bridge/local3_vs_hybrid3_ranking.md` |
| hybrid3 family evidence sensitivity | `results/prl_reference_bridge/family_evidence_sensitivity.md` |
| bridge gamma projection | `results/prl_reference_bridge/gamma_projection.md` |
| bridge frequency-cut diagnostic | `results/prl_reference_bridge/frequency_cut_evidence.md` |
| bridge robustness and QMC/TI checks | `results/prl_reference_bridge/robustness_budget.csv`, `results/prl_reference_bridge/ti_crosscheck_top_models.csv` |

## Minimal Verification Commands

Run from the project root:

```bash
python3 code/gwb_templates.py
python3 -m py_compile \
  code/bayes_factors_ptarcade_groundtruth.py \
  code/bayes_factors_ptarcade_ceffyl_density.py \
  code/bayes_factors_ptarcade_ceffyl_density_extended.py \
  code/bayes_factors_Nbins_scan.py \
  code/prl_official_density_robustness.py \
  code/prl_car_null_calibration.py \
  code/prl_H4_H5_H6_experiments.py \
  code/prl_H4_env_robustness.py \
  code/prl_H7_astrophysical_family.py \
  code/prl_decisive_evidence_figure.py \
  code/prl_bridge_evidence_figure.py \
  code/prl_reference_bridge_pipeline.py \
  code/prl_H9_bin_driver_analysis.py \
  code/prl_H10_systematic_envelope.py \
  code/prl_evidence_ti_qmc_crosscheck.py \
  code/prl_package_static_gate.py
```

The template sanity-check command should finish with `ALL CHECKS PASSED`.

To reproduce the fast official-density evidence checks:

```bash
python3 code/bayes_factors_ptarcade_ceffyl_density.py
python3 code/bayes_factors_ptarcade_ceffyl_density_extended.py
```

The expected PTArcade-BHB-prior-baseline log-Bayes factors are:

| Model | Expected lnB |
|---|---:|
| SIGW-Gauss | `+4.225 ± 0.173` |
| SIGW-delta | `+3.647 ± 0.169` |
| Cosmic superstrings (`super.py`) | `+3.338 ± 0.168` |

Do not include `results/T2_NG15yr/chains_production_HD_hot/` in convergence
diagnostics. That chain was diagnosed as invalid because it remained at
`[-inf, -inf, 0, 1]` metadata with zero accepted jumps.

To run the PRL package static gate:

```bash
python3 code/prl_evidence_ti_qmc_crosscheck.py
python3 code/prl_decisive_evidence_figure.py
python3 code/prl_reference_bridge_pipeline.py all --profile production --force
python3 code/prl_bridge_evidence_figure.py
python3 code/prl_H10_systematic_envelope.py
(cd theory && tectonic --keep-logs --keep-intermediates --outdir pdf/revtex paper_prl_submission.tex)
(cd theory && tectonic --keep-logs --keep-intermediates --outdir pdf/revtex prl_supplement.tex)
python3 code/prl_package_static_gate.py
```

As of 2026-04-23, Tectonic is available locally and the REVTeX main text and
supplement compile to `theory/pdf/revtex/`.  The public code/data release is
archived at Zenodo DOI `10.5281/zenodo.19688587` and mirrored at
`https://github.com/duan0004/paperdata`, release `v1.0.1`, commit
`2bc1dfeee7d452f6c3b58a8fd8f339e6ff0091f6`.  The cross-PTA bridge update is
included in the repository history after that release and should receive its
own Zenodo snapshot before final journal submission.
