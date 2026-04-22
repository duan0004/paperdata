# T1b.1 Pilot Report — Bin-to-bin posterior covariance of NG15yr HD 14-bin free spectrum

*Generated:* 2026-04-21 08:55:32
*Source:* `data/NG15yr/tutorials/presampled_cores/hd_30f_fs.core`
*Post-burn samples:* 9,188

## Pilot gates (kill-switch 2026-05-05)

- **G1 (structure)**: ≥3 signal modes above Marchenko-Pastur noise AND effective rank ≥3 → **FAIL**
    - #signal modes (λ > 1.2 × 1.0796) = **1**
    - effective rank r_eff = (Σλ)²/Σλ² = **13.56**
    - #|corr_ij|≥0.1 off-diagonal pairs = **7**, max |off-diagonal| = **0.227**

- **G2 (Gaussianity)**: ≥12/14 bins pass |skew|<0.5 & |kurt-3|<2 → **FAIL**
    - 8/14 bins pass
    - skew range [-26.237, +0.042]
    - kurt range [1.691, 1252.546]

## Verdict: **INVESTIGATE — see failure modes below**

## Top-5 eigenvalues of the correlation matrix

| k | λ_k | participation ratio |
|---|-----|---------------------|
| 1 | 1.4826 | 4.49 |
| 2 | 1.2133 | 7.84 |
| 3 | 1.0843 | 4.50 |
| 4 | 1.0796 | 5.86 |
| 5 | 1.0172 | 3.20 |

## Interpretation

The correlation matrix is normalized so a diagonal ceffyl refit would see eigenvalues all ≡ 1.
Eigenvalues > 1 indicate *positively correlated* bin groups (a coherent inflation of several bins);
eigenvalues < 1 indicate anti-correlated directions — these are the directions a wrong-shape template
has to fit against, and are the **discriminating** modes for CAR.  Ceffyl effectively inflates the
error along small-λ directions (treating them as σ²=1 instead of σ²=λ<1), wasting constraining power.

## Files

- `mu_bin.npy`, `sigma_bin.npy`, `corr_bin.npy` — the 14-D summary statistics
- `eigvals.npy`, `eigvecs.npy` — correlation-matrix spectrum
- `eigenspectrum.png` — figure (a) spectrum, (b) correlation heatmap, (c) top-3 eigenvectors
- `pilot_report.json` — machine-readable dump