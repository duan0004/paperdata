#!/usr/bin/env python3
"""
extract_bin_covariance.py — Plan D / Phase 1b / T1b.1
======================================================

Extract the 14×14 bin-to-bin posterior covariance Σ_bin of {log10_ρ_k}
from the public NANOGrav 15 yr HD free-spectrum core (`hd_30f_fs.core`),
using only the 14 lowest Fourier bins per the NG15yr-NewPhysics §III.B
convention adopted throughout this project.

Motivation (Plan D §I):
  Standard ceffyl-style refits treat the 14-bin free-spectrum posterior
  as a product of 14 independent per-bin KDEs.  §IV.C of the paper
  diagnoses that the discarded off-diagonal covariance is plausibly the
  *dominant* source of the 2-dex absolute-Bayes-factor gap w.r.t.
  NG15yr Table 3.

  Covariance-Aware Refit (CAR) inverts this: use the true 14-D joint
  posterior (via its mean μ_bin and covariance Σ_bin, adequate if the
  marginal is approximately Gaussian in log10_ρ; verified below) as the
  likelihood kernel, and project the spectral templates onto the
  eigenmodes of Σ_bin.  Eigenmodes with small eigenvalues are the
  discriminating directions — they cost a lot of lnL for a model whose
  spectrum doesn't match the data along that axis.

Pilot acceptance (kill-switch gate 2026-05-05):
  G1. At least 3 eigenvalues of the correlation matrix
      (normalized so diagonal = 1) carry non-trivial structure,
      i.e. |corr_ij| ≥ 0.1 for at least 3 independent (i,j) pairs, OR
      equivalently the effective rank r_eff = (Σλ)² / Σλ² is ≥ 3.
  G2. Joint posterior is approximately Gaussian in log10_ρ:
      per-bin |skewness| < 0.5 and |kurtosis - 3| < 2 for at least 12/14
      bins.  (If violated, CAR needs a full-KDE multivariate kernel,
      not Gaussian — still feasible, flag for T1b.3.)

Outputs:
  results/T2_NG15yr/covariance/mu_bin.npy           14-vector
  results/T2_NG15yr/covariance/sigma_bin.npy        14×14
  results/T2_NG15yr/covariance/corr_bin.npy         14×14 (Pearson)
  results/T2_NG15yr/covariance/eigvals.npy          14-vector (desc)
  results/T2_NG15yr/covariance/eigvecs.npy          14×14 (columns)
  results/T2_NG15yr/covariance/pilot_report.json    numerical summary
  results/T2_NG15yr/covariance/pilot_report.md      human-readable
  results/T2_NG15yr/covariance/eigenspectrum.png    figure

Data provenance:
  arXiv:2306.16213 tutorial data — `hd_30f_fs.core`
"""
import os, sys, json, time
import numpy as np
import h5py
from scipy import stats

# ── paths ───────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORE = os.path.join(BASE, "data", "NG15yr", "tutorials", "presampled_cores", "hd_30f_fs.core")
OUT  = os.path.join(BASE, "results", "T2_NG15yr", "covariance")
os.makedirs(OUT, exist_ok=True)

N_BINS     = 14             # NG15yr-NewPhysics §III.B convention
T_OBS_YR   = 16.03          # NG15yr observing baseline
T_OBS_S    = T_OBS_YR * 3.15576e7

print("=" * 70)
print("T1b.1  Extract 14×14 bin-to-bin posterior covariance")
print("=" * 70)

# ── load hd_30f_fs.core ─────────────────────────────────────────────
t0 = time.time()
with h5py.File(CORE, "r") as f:
    params_raw = f["params"][...]
    params = [p.decode() if isinstance(p, bytes) else p for p in params_raw]
    burn = int(f["metadata/burn"][()])
    chain_full = f["chain"][burn:, :]
print(f"  loaded core: {chain_full.shape[0]:,} post-burn samples × {chain_full.shape[1]} params ({time.time()-t0:.1f}s)")

# locate the 14 log10_ρ_k columns
idx = []
for k in range(N_BINS):
    name = f"gw_hd_log10_rho_{k}"
    if name not in params:
        print(f"  ERROR: column {name!r} not found in core"); sys.exit(1)
    idx.append(params.index(name))
rho = chain_full[:, idx]            # shape (N_samples, 14)
print(f"  extracted rho block: {rho.shape}")

# ── Gaussianity check (G2) ──────────────────────────────────────────
skew = stats.skew(rho, axis=0)
kurt = stats.kurtosis(rho, axis=0, fisher=False)   # NOT excess kurtosis
n_gauss = int(np.sum((np.abs(skew) < 0.5) & (np.abs(kurt - 3.0) < 2.0)))
print(f"\n  G2 Gaussianity: {n_gauss}/{N_BINS} bins pass |skew|<0.5 & |kurt-3|<2")
print(f"       skew range [{skew.min():+.3f}, {skew.max():+.3f}]")
print(f"       kurt range [{kurt.min():.3f}, {kurt.max():.3f}]")

# ── covariance + correlation ────────────────────────────────────────
mu    = rho.mean(axis=0)
sigma = np.cov(rho, rowvar=False)
d     = np.sqrt(np.diag(sigma))
corr  = sigma / np.outer(d, d)

# zero-out the diagonal for off-diagonal-mass accounting
off = corr - np.eye(N_BINS)
n_strong = int(np.sum(np.abs(off) >= 0.1) // 2)   # symmetric → halve
max_offdiag = float(np.max(np.abs(off)))

# ── eigendecomposition of the correlation matrix ────────────────────
# (correlation, not covariance — removes per-bin scale so the eigenvalues
#  tell us about *shape* correlations, not overall amplitude)
eigvals, eigvecs = np.linalg.eigh(corr)
order   = np.argsort(eigvals)[::-1]        # descending
eigvals = eigvals[order]
eigvecs = eigvecs[:, order]

# effective rank = (Σλ)² / Σλ²
r_eff = float((eigvals.sum())**2 / np.sum(eigvals**2))

# participation ratio by mode — how many bins contribute per mode
part_ratio = 1.0 / np.sum(eigvecs**4, axis=0)

# ── G1 pilot criterion ──────────────────────────────────────────────
# "At least 3 eigenvalues carry non-trivial structure":
# declare top-k modes nontrivial if their eigenvalues are ≥ 1.2× the
# noise-floor expectation for a 14-dim PCA on N samples.  The
# Marchenko-Pastur upper bound for an IID-null is λ_+ = (1 + √(p/N))²
# with p=14, N=len(rho).
N_samp = rho.shape[0]
p = N_BINS
lam_noise_upper = (1.0 + np.sqrt(p / N_samp))**2
n_signal_modes = int(np.sum(eigvals > 1.2 * lam_noise_upper))

g1_pass = (n_signal_modes >= 3) and (r_eff >= 3.0)
g2_pass = (n_gauss >= 12)

print(f"\n  eigenvalues (desc): {np.array2string(eigvals, precision=3, separator=', ')}")
print(f"  top-5 eigvals:       {eigvals[:5]}")
print(f"  MP noise upper λ_+:  {lam_noise_upper:.4f}  (1.2× = {1.2*lam_noise_upper:.4f})")
print(f"  #signal modes (λ > 1.2 λ_+): {n_signal_modes}")
print(f"  effective rank r_eff:        {r_eff:.2f}")
print(f"  #|corr_ij|≥0.1 pairs (i<j):  {n_strong}")
print(f"  max |off-diagonal corr|:     {max_offdiag:.3f}")

print(f"\n  ── Pilot gates ──")
print(f"    G1 (≥3 signal modes AND r_eff ≥ 3): {'PASS' if g1_pass else 'FAIL'}")
print(f"    G2 (≥12/14 bins approx Gaussian):   {'PASS' if g2_pass else 'FAIL'}")

# ── save numerical outputs ──────────────────────────────────────────
np.save(os.path.join(OUT, "mu_bin.npy"),    mu)
np.save(os.path.join(OUT, "sigma_bin.npy"), sigma)
np.save(os.path.join(OUT, "corr_bin.npy"),  corr)
np.save(os.path.join(OUT, "eigvals.npy"),   eigvals)
np.save(os.path.join(OUT, "eigvecs.npy"),   eigvecs)

report = {
    "task": "T1b.1",
    "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
    "source_core": os.path.relpath(CORE, BASE),
    "N_samples_postburn": int(N_samp),
    "N_bins": N_BINS,
    "T_obs_yr": T_OBS_YR,
    "mu_bin": mu.tolist(),
    "diag_sigma": np.diag(sigma).tolist(),
    "gaussianity": {
        "skew": skew.tolist(),
        "kurtosis": kurt.tolist(),
        "n_bins_pass": n_gauss,
    },
    "correlation_structure": {
        "n_strong_offdiag_pairs_|corr|>=0.1": n_strong,
        "max_abs_offdiag": max_offdiag,
    },
    "eigendecomposition": {
        "eigvals_desc": eigvals.tolist(),
        "effective_rank": r_eff,
        "participation_ratio": part_ratio.tolist(),
        "marchenko_pastur_noise_upper": lam_noise_upper,
        "n_signal_modes_above_1p2_MP": n_signal_modes,
    },
    "pilot_gates": {
        "G1_pass": bool(g1_pass),
        "G2_pass": bool(g2_pass),
        "go_no_go": "GO" if (g1_pass and g2_pass) else "INVESTIGATE",
    },
}
with open(os.path.join(OUT, "pilot_report.json"), "w") as f:
    json.dump(report, f, indent=2)

# ── figure ──────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# (a) eigenvalue spectrum + MP noise floor
ax = axes[0]
ax.plot(range(1, N_BINS+1), eigvals, "o-", color="C0", label="Σ_bin correlation eigenvalues")
ax.axhline(lam_noise_upper, ls="--", color="grey", label=f"MP upper λ_+={lam_noise_upper:.3f}")
ax.axhline(1.2*lam_noise_upper, ls=":", color="k", label="1.2× noise (signal cut)")
ax.set_xlabel("mode index (desc)"); ax.set_ylabel("eigenvalue")
ax.set_yscale("log")
ax.set_title(f"(a) Eigenspectrum — {n_signal_modes} signal modes, r_eff={r_eff:.2f}")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# (b) correlation heatmap
ax = axes[1]
im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, origin="lower")
ax.set_title(f"(b) Bin-to-bin correlation matrix (max |off-diag|={max_offdiag:.2f})")
ax.set_xlabel("bin $k$"); ax.set_ylabel("bin $k'$")
plt.colorbar(im, ax=ax, shrink=0.85)

# (c) top-3 eigenvectors
ax = axes[2]
for k in range(min(3, N_BINS)):
    ax.plot(range(N_BINS), eigvecs[:, k], "o-", label=f"mode {k+1} (λ={eigvals[k]:.3f})")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("bin index"); ax.set_ylabel("eigenvector coefficient")
ax.set_title("(c) Top-3 eigenmodes (structure of the discarded covariance)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(OUT, "eigenspectrum.png"), dpi=130)
print(f"\n  figure saved: {os.path.relpath(os.path.join(OUT, 'eigenspectrum.png'), BASE)}")

# ── human-readable report ───────────────────────────────────────────
md = [
    "# T1b.1 Pilot Report — Bin-to-bin posterior covariance of NG15yr HD 14-bin free spectrum",
    "",
    f"*Generated:* {time.strftime('%Y-%m-%d %H:%M:%S')}",
    f"*Source:* `{os.path.relpath(CORE, BASE)}`",
    f"*Post-burn samples:* {N_samp:,}",
    "",
    "## Pilot gates (kill-switch 2026-05-05)",
    "",
    f"- **G1 (structure)**: ≥3 signal modes above Marchenko-Pastur noise AND effective rank ≥3 → **{'PASS' if g1_pass else 'FAIL'}**",
    f"    - #signal modes (λ > 1.2 × {lam_noise_upper:.4f}) = **{n_signal_modes}**",
    f"    - effective rank r_eff = (Σλ)²/Σλ² = **{r_eff:.2f}**",
    f"    - #|corr_ij|≥0.1 off-diagonal pairs = **{n_strong}**, max |off-diagonal| = **{max_offdiag:.3f}**",
    "",
    f"- **G2 (Gaussianity)**: ≥12/14 bins pass |skew|<0.5 & |kurt-3|<2 → **{'PASS' if g2_pass else 'FAIL'}**",
    f"    - {n_gauss}/{N_BINS} bins pass",
    f"    - skew range [{skew.min():+.3f}, {skew.max():+.3f}]",
    f"    - kurt range [{kurt.min():.3f}, {kurt.max():.3f}]",
    "",
    f"## Verdict: **{'GO — proceed with Plan D' if (g1_pass and g2_pass) else 'INVESTIGATE — see failure modes below'}**",
    "",
    "## Top-5 eigenvalues of the correlation matrix",
    "",
    "| k | λ_k | participation ratio |",
    "|---|-----|---------------------|",
]
for k in range(5):
    md.append(f"| {k+1} | {eigvals[k]:.4f} | {part_ratio[k]:.2f} |")
md.append("")
md.append("## Interpretation")
md.append("")
md.append("The correlation matrix is normalized so a diagonal ceffyl refit would see eigenvalues all ≡ 1.")
md.append("Eigenvalues > 1 indicate *positively correlated* bin groups (a coherent inflation of several bins);")
md.append("eigenvalues < 1 indicate anti-correlated directions — these are the directions a wrong-shape template")
md.append("has to fit against, and are the **discriminating** modes for CAR.  Ceffyl effectively inflates the")
md.append("error along small-λ directions (treating them as σ²=1 instead of σ²=λ<1), wasting constraining power.")
md.append("")
md.append("## Files")
md.append("")
md.append("- `mu_bin.npy`, `sigma_bin.npy`, `corr_bin.npy` — the 14-D summary statistics")
md.append("- `eigvals.npy`, `eigvecs.npy` — correlation-matrix spectrum")
md.append("- `eigenspectrum.png` — figure (a) spectrum, (b) correlation heatmap, (c) top-3 eigenvectors")
md.append("- `pilot_report.json` — machine-readable dump")

with open(os.path.join(OUT, "pilot_report.md"), "w") as f:
    f.write("\n".join(md))

print(f"\n  report saved: {os.path.relpath(os.path.join(OUT, 'pilot_report.md'), BASE)}")
print(f"\n[T={time.time()-t0:.1f}s] done.")
print(f"\n>>> Pilot verdict: {'GO' if (g1_pass and g2_pass) else 'INVESTIGATE'} <<<")
