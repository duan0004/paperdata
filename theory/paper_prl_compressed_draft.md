# Curved SMBHB Spectra Compete with New-Physics Explanations in PTA Source Identification

**Ran DUAN**  
**Affiliation**: National Astronomical Observatories, Chinese Academy of Sciences  
**Email**: duanran@nao.cas.cn

**Code/data release**: Zenodo DOI `10.5281/zenodo.19688472`; GitHub mirror https://github.com/duan0004/paperdata, release `v1.0.0`, commit `7a602df4d0369dea09b3a18d526a1468bbc436d0`.

## Abstract

The NANOGrav 15-yr Hellings--Downs signal admits both astrophysical and cosmological interpretations. We test whether the published new-physics Bayes factors remain discriminating once template provenance, density construction, and curved-SMBHB controls are matched. Official PTArcade templates and official PTArcade `ceffyl` densities reproduce the published log-Bayes factors for SIGW-Gaussian, SIGW-delta, and cosmic superstrings at sub-nat agreement. On that calibrated scale, one environmental SMBHB model and two phenomenological low-frequency curvature surrogates yield $\ln\mathcal{B}=3.45$--$3.84$, placing three distinct curved-SMBHB parameterizations within one nat of the leading cosmological evidence tier. The $0.681$-nat gap between the leading SIGW-Gaussian model and the best curved-SMBHB control is comparable to the identified calibration and prior-boundary systematics. Current data therefore do not yet support robust source discrimination between leading cosmological stochastic spectra and calibrated curved-SMBHB controls.

## I. Introduction

Pulsar timing arrays have opened the nanohertz gravitational-wave window. The NANOGrav 15-year data set, with 67 millisecond pulsars over $T=16.03\,\mathrm{yr}$, reports a common red process with Hellings--Downs spatial correlations \cite{NANOGrav_15yr_GWB}. The companion noise analysis \cite{NANOGrav_15yr_noise} provides the white-noise budget used in the public analyses. The same signal is compatible with an astrophysical supermassive-black-hole-binary (SMBHB) background and with several cosmological sources.

The astrophysical baseline is a circular, GW-driven SMBHB population with characteristic strain $h_c=A(f/f_\mathrm{ref})^{-2/3}$, corresponding to the PTA power-spectral index $\gamma=13/3$ \cite{Phinney_2001}. NANOGrav finds $A_{\mathrm{HD},13/3}=2.4^{+0.7}_{-0.6}\times10^{-15}$ for fixed $\gamma=13/3$, while the free-index fit gives $\gamma_\mathrm{HD}\simeq3.2\pm0.6$ \cite{NANOGrav_15yr_GWB}. The new-physics search \cite{NANOGrav_15yr_NewPhysics}, using PTArcade \cite{PTArcade_2023}, reports strong-but-not-decisive Bayes factors relative to the SMBHB baseline for several leading stochastic alternatives, including SIGW-Gaussian ($\mathcal{B}=57$), SIGW-delta ($\mathcal{B}=44$), and cosmic superstrings ($\mathcal{B}=46$), under the Kass--Raftery convention \cite{Kass_Raftery_1995}. We show below that one environmental SMBHB model and two phenomenological curved-SMBHB surrogates are competitive with these alternatives on the same evidence scale.

The central problem is reproducibility. Order-unity shifts in $\ln\mathcal{B}$ can arise from the scalar-induced-GW kernel, the $u+v=\sqrt{3}$ radiation-era resonance \cite{Kohri_Terada_2018}, the prior boundary on the peak frequency, or the density estimate used for the free-spectrum likelihood. We therefore separate three questions: whether our PTA pipeline reproduces the published HD posterior, whether our local spectral library reproduces the qualitative ranking, and whether official PTArcade spectra plus official density products reproduce the published absolute Bayes factors.

## II. Templates and PTA Likelihood

We model an isotropic stochastic background through the one-sided strain power spectral density $S_h(f)$ and
$$
\Omega_\mathrm{GW}(f)=\frac{2\pi^2}{3H_0^2}f^3S_h(f),
$$
with $H_0=67.4\,\mathrm{km\,s^{-1}\,Mpc^{-1}}$. The template library includes: SMBHB power law; SMBHB with environmental hardening \cite{Sesana_2013,Sampson_2015}; SIGW-Gaussian; SIGW-delta; stable Nambu--Goto strings; cosmic superstrings; and two first-order phase-transition spectra, bubble collision and sound waves. All $\Omega_\mathrm{GW}$ calculations are dimension-checked in `code/gwb_templates.py`.

For the PTA likelihood we use ENTERPRISE \cite{ENTERPRISE_2019}. White-noise parameters EFAC, EQUAD, and ECORR are fixed to the official NANOGrav 15-yr noise dictionary \cite{NANOGrav_15yr_noise}. Per-pulsar red noise is modeled as
$$
P_a(f)=\frac{A_a^2}{12\pi^2}\left(\frac{f}{f_\mathrm{ref}}\right)^{-\gamma_a}f_\mathrm{ref}^{-3},
\qquad
\phi_{a,i}=P_a(f_i)/T,
$$
with 30 red-noise Fourier components per pulsar. The common process uses 14 GWB Fourier bins at $f_i=i/T$, matching the NANOGrav low-frequency truncation \cite{NANOGrav_15yr_GWB,NANOGrav_15yr_NewPhysics}. Timing-model offsets are marginalized with `MarginalizingTimingModel`; `.par` files are not modified.

The full HD production model contains 140 sampled parameters: 67 intrinsic-red-noise amplitudes, 67 intrinsic-red-noise spectral indices, and the two common-process parameters $(\log_{10}A_\mathrm{GWB},\gamma_\mathrm{GWB})$, with HD spatial correlations active from the first step. Model evidences are computed either from local free-spectrum KDE products or from the official PTArcade `ceffyl` density grid using `dynesty` nested sampling with a target `dlogz=0.1`.

## III. Pipeline Validation

We first validate against the public NANOGrav posterior products. Direct ingestion of the pre-sampled HD chain gives $\log_{10}A_\mathrm{HD,vg}=-14.20\pm0.13$ and $\gamma=3.25\pm0.35$, matching the published result \cite{NANOGrav_15yr_GWB}. This confirms the data loading, noise dictionary, and posterior parsing.

We then run our own HD PTMCMC chain. A $2\times10^6$-step warm chain gives
$$
\log_{10}A_\mathrm{GWB}=-14.199\pm0.122,\qquad
\gamma_\mathrm{GWB}=3.272\pm0.321.
$$
A second $5\times10^5$-step chain initialized from dispersed cold conditions, $(\log_{10}A,\gamma)=(-15.5,5.0)$, converges to $\log_{10}A=-14.251\pm0.168$ and $\gamma=3.401\pm0.418$. Tail-aligned warm+cold segments yield $\hat R(\log_{10}A)=1.047$ and $\hat R(\gamma)=1.047$, below an a priori convergence threshold of $1.05$. A third hot-chain attempt was rejected because it remained at $-\infty$ posterior metadata with zero accepted jumps.

## IV. Bayes-Factor Checks

Our local `ceffyl`-style refit starts from the released HD 30-bin free-spectrum core, builds independent per-bin KDEs for the 14 lowest bins, and evaluates each spectral template through
$$
\ln\mathcal{L}(\theta)=\sum_{k=1}^{14}
\ln\mathrm{KDE}_k\!\left[\log_{10}\rho_k^\mathrm{model}(\theta)\right].
$$
This local refit reproduces the qualitative NANOGrav ranking but not the published absolute values. With our semi-analytic templates, SIGW-Gaussian gives $\ln\mathcal{B}=-1.59\pm0.15$ relative to the fixed-$\gamma$ SMBHB baseline, while the deliberately flat stable-string proxy is strongly disfavored. Sensitivity sweeps show that KDE bandwidth, not nested-sampler seed or `nlive`, dominates this local-refit uncertainty.

We next isolate the source of the discrepancy. The released HD free-spectrum covariance is nearly diagonal: max $|\mathrm{corr}_{ij}|=0.227$, mean $|\mathrm{corr}_{ij}|=0.030$, and effective rank $13.56/14$. A CAR projection gives a negligible correction for SIGW-Gaussian versus SMBHB, $\Delta\ln\mathcal{Z}\simeq-0.006$ nat, and a modest correction for the poor flat-string proxy, $\simeq+1.31$ nat. Covariance restoration therefore cannot explain a multi-nat absolute-Bayes-factor gap.

The decisive test uses the official NANOGrav/PTArcade model files and the official PTArcade `ceffyl` density product. PTArcade model files return $h^2\Omega_\mathrm{GW}$, so the residual conversion is
$$
\rho^2=\frac{(H_0/h)^2\,h^2\Omega_\mathrm{GW}}
{8\pi^4 f^5 T}.
$$
This conversion was checked against `ptarcade.models_utils.omega2cross` to relative precision below $5\times10^{-16}$.

The official-density robustness results, averaged over five seed/live-point configurations, are:

| Model | This work mean $\ln\mathcal{B}$ | Published $\ln\mathcal{B}$ | Mean residual |
|---|---:|---:|---:|
| SIGW-Gaussian | $+4.520\pm0.185$ | $\ln(57)=+4.043$ | $+0.477$ nat |
| SIGW-delta | $+3.930\pm0.199$ | $\ln(44)=+3.784$ | $+0.146$ nat |
| Cosmic superstrings (`super.py`) | $+3.558\pm0.145$ | $\ln(46)=+3.829$ | $-0.271$ nat |
| Eccentricity-inspired curved SMBHB | $+3.839\pm0.043$ | not reported | best tested representative |
| Broken-power-law curved SMBHB | $+3.656\pm0.154$ | not reported | best tested representative |
| Environmental SMBHB | $+3.453\pm0.062$ | not reported | best tested representative |

For the three published cosmological templates, the largest absolute residual across the five configurations is $0.694$ nat. In contrast, replacing only the SIGW-Gaussian spectrum while retaining our local scipy-KDE density gives $\ln\mathcal{B}=+0.687\pm0.147$, still low by $3.36$ nats. The absolute evidence alignment is therefore controlled primarily by official density construction and template fidelity, not by bin covariance. The curved-SMBHB controls show that the result is not a single-model accident: the best curved-SMBHB model lies only $0.681$ nat below SIGW-Gaussian and $0.091$ nat below SIGW-delta. The $0.681$-nat gap is comparable to the identified calibration and prior-boundary systematics, so it should not be interpreted as robust source discrimination. Independent prior-space Sobol thermodynamic integration agrees with the dynesty robust means to within $0.12$ nat for the leading cosmological templates and within $0.063$ nat for the tested curved-SMBHB representatives.

The Supplemental Material gives the full official-density stochastic-model ranking and prior-boundary tests. In that scan the leading six models are SIGW-Gaussian, SIGW-delta, SIGW-box, cosmic superstrings, metastable strings, and domain walls. Extending the SIGW-delta upper prior from $\log_{10}f_\mathrm{peak}\le -5$ to $-3$ increases its evidence by $0.426$ nat, a modest but non-negligible boundary systematic.

As an orthogonal diagnostic, we scan the number of retained Fourier bins. The flat stable-string proxy drifts monotonically from $\Delta\ln\mathcal{B}=-24.25$ to $-28.12$ to $-29.16$ for $N_\mathrm{bins}=6,14,30$, whereas the SIGW-Gaussian proxy does not show a strong monotonic trend. We retain $\ln\mathcal{B}(N_\mathrm{bins})$ as a spectral-mismatch diagnostic, not as a physical exclusion of strings.

## V. Interpretation

The NANOGrav 15-yr data do not yet identify the physical origin of the background. On a calibrated evidence scale, current leading new-physics preferences are compressed into the same model-comparison tier as curved SMBHB spectra. The most robust near-term discriminant is spectral shape at low frequency: SMBHBs give $h^2\Omega_\mathrm{GW}\propto f^{2/3}$, whereas causal cosmological sources such as SIGW and first-order phase transitions rise as $f^3$ below their characteristic scale. The cumulative-bin scan makes this quantitative: the lowest eight to ten Fourier bins already determine the matched-scale ranking to within $0.5$ nat for the leading stochastic templates.

The practical lesson is methodological. A local free-spectrum KDE product is useful for fast diagnostics, but it is not automatically equivalent to the official PTArcade `ceffyl` density. For absolute Bayes factors at the sub-nat level, the official density product, the exact template implementation, and the matched SMBHB prior baseline must be used together; cosmological alternatives must also be calibrated against curved SMBHB controls. CAR/bin covariance is a useful diagnostic but is not the dominant error budget for the leading SIGW/string comparisons.

## VI. Conclusions

We constructed and validated a spectral-template analysis pipeline for the NANOGrav 15-year background. The ENTERPRISE/PTMCMCSampler HD run reproduces the published common-process posterior and passes a two-chain $\hat R<1.05$ diagnostic. Local refits reproduce the qualitative model ranking but expose sensitivity to KDE density construction. Official PTArcade templates plus official PTArcade `ceffyl` densities reproduce the published SIGW-Gaussian, SIGW-delta, and cosmic-superstring log-Bayes factors at sub-nat precision across seed/live-point repetitions, while a curved SMBHB family overlaps the leading cosmological evidence range. The PRL scientific gate is therefore stronger than a calibration note: current data do not robustly discriminate leading cosmological stochastic spectra from calibrated curved-SMBHB controls. The code/data release is archived at Zenodo DOI `10.5281/zenodo.19688472`; the remaining work before PRL submission is final human approval of any acknowledgments.

## References

Use the reference set in `theory/paper_full_draft.md`; priority citations are
\cite{NANOGrav_15yr_GWB,NANOGrav_15yr_noise,NANOGrav_15yr_NewPhysics,PTArcade_2023,Phinney_2001,Sesana_2013,Sampson_2015,Kass_Raftery_1995,ENTERPRISE_2019,PTMCMCSampler_2017,Speagle_2020,Kohri_Terada_2018,Caprini_Figueroa_2018}.
