#!/usr/bin/env python3
"""
extract_bin_covariance_lowk.py — Plan D diagnostic follow-up
==============================================================

After the full-14-bin pilot (extract_bin_covariance.py) reported
max |off-diag corr| = 0.227 and only 1 signal eigenmode, we investigate
whether the correlation structure concentrates in the signal-dominated
low-k subspace.

Rationale:
  NG15yr HD signal is concentrated in bins 0-7 (~1-15 nHz).  High-k bins
  (8-13) are noise-dominated, their marginals are rail-truncated against
  the log10_ρ prior floor, and they cannot correlate strongly because
  each is consistent with zero independently.  The bin-independence
  approximation is *automatically correct* in the noise-dominated
  subspace — so the correlation signal, if any, must live in low-k.

Test:
  Rerun the 14-bin pilot restricted to the first N_sig ∈ {6, 8, 10} bins.
  Report eigenspectrum + correlation strength in each subspace.

Output:
  results/T2_NG15yr/covariance/lowk_investigation.json
  results/T2_NG15yr/covariance/lowk_eigenspectrum.png
"""
import os, sys, json, time
import numpy as np
import h5py
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORE = os.path.join(BASE, "data", "NG15yr", "tutorials", "presampled_cores", "hd_30f_fs.core")
OUT  = os.path.join(BASE, "results", "T2_NG15yr", "covariance")

N_TOTAL = 14
N_SIG_LIST = [4, 6, 8, 10, 12, 14]

with h5py.File(CORE, "r") as f:
    params = [p.decode() if isinstance(p, bytes) else p for p in f["params"][...]]
    burn = int(f["metadata/burn"][()])
    chain = f["chain"][burn:, :]

idx = [params.index(f"gw_hd_log10_rho_{k}") for k in range(N_TOTAL)]
rho = chain[:, idx]
N_samp = rho.shape[0]
print(f"loaded rho block: {rho.shape}\n")

results = []
for n in N_SIG_LIST:
    sub = rho[:, :n]
    skew = stats.skew(sub, axis=0)
    kurt = stats.kurtosis(sub, axis=0, fisher=False)
    n_gauss = int(np.sum((np.abs(skew) < 0.5) & (np.abs(kurt - 3.0) < 2.0)))

    sigma = np.cov(sub, rowvar=False)
    d = np.sqrt(np.diag(sigma))
    corr = sigma / np.outer(d, d)
    off = corr - np.eye(n)
    n_strong = int(np.sum(np.abs(off) >= 0.1) // 2)
    max_off  = float(np.max(np.abs(off)))
    mean_abs_off = float(np.mean(np.abs(off[np.triu_indices(n, k=1)])))

    eigvals = np.sort(np.linalg.eigvalsh(corr))[::-1]
    lam_mp = (1.0 + np.sqrt(n / N_samp))**2
    n_sig_modes = int(np.sum(eigvals > 1.2 * lam_mp))
    r_eff = float(eigvals.sum()**2 / np.sum(eigvals**2))
    # anisotropy: ratio of largest to smallest eigenvalue
    aniso = float(eigvals[0] / eigvals[-1])

    row = {
        "n_bins_used": n,
        "eigvals_top5": eigvals[:min(5, n)].tolist(),
        "eigval_max": float(eigvals[0]),
        "eigval_min": float(eigvals[-1]),
        "anisotropy_max_over_min": aniso,
        "MP_noise_upper": lam_mp,
        "n_signal_modes_above_1p2_MP": n_sig_modes,
        "r_eff_rel_to_n": r_eff / n,   # 1.0 = isotropic, <1 = structured
        "n_strong_offdiag_pairs_abs_corr_ge_0p1": n_strong,
        "max_abs_offdiag_corr": max_off,
        "mean_abs_offdiag_corr": mean_abs_off,
        "n_bins_gaussian": n_gauss,
    }
    results.append(row)
    print(f"  n={n:2d} | #sig modes={n_sig_modes} | max|corr|={max_off:.3f}"
          f" | mean|corr|={mean_abs_off:.3f}"
          f" | λ_max={eigvals[0]:.3f} | aniso={aniso:.2f}"
          f" | gauss={n_gauss}/{n}")

# ── decision logic ─────────────────────────────────────────────────
# Recast G1: plan D needs EITHER (i) ≥3 signal modes on full 14 bins,
# OR (ii) ≥2 signal modes on the low-k signal subspace n∈{6,8,10} with
# max|corr| ≥ 0.3 (i.e. stronger local structure in the signal bins).
lowk_sig_cases = [r for r in results if r["n_bins_used"] in (6, 8, 10)]
plan_d_viable_by_lowk = any(
    (r["n_signal_modes_above_1p2_MP"] >= 2) and (r["max_abs_offdiag_corr"] >= 0.3)
    for r in lowk_sig_cases
)
plan_d_viable_by_full = any(
    r["n_bins_used"] == 14 and r["n_signal_modes_above_1p2_MP"] >= 3
    for r in results
)

# Softer criterion: even without classic "signal modes", if λ_max ≥ 1.3
# on the low-k subspace, the CAR has enough structure to provide ~1-nat
# lnB corrections for narrow-feature templates — this keeps Plan D viable
# but at "modest" rather than "dominant" magnitude.
soft_viable = any(
    r["n_bins_used"] in (6, 8, 10) and r["eigval_max"] >= 1.3
    for r in results
)

verdict = {
    "plan_d_viable_by_lowk_hard":  bool(plan_d_viable_by_lowk),
    "plan_d_viable_by_full_14":    bool(plan_d_viable_by_full),
    "plan_d_viable_soft_lambda_max_ge_1p3": bool(soft_viable),
}

print("\n── Revised viability assessment ──")
print(f"  Hard (low-k ≥2 signal modes AND max|corr|≥0.3): "
      f"{'PASS' if plan_d_viable_by_lowk else 'FAIL'}")
print(f"  Full 14-bin (≥3 signal modes):                   "
      f"{'PASS' if plan_d_viable_by_full else 'FAIL'}")
print(f"  Soft (low-k λ_max ≥ 1.3):                         "
      f"{'PASS' if soft_viable else 'FAIL'}")

out = {
    "task": "T1b.1-followup",
    "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
    "scan": results,
    "verdict": verdict,
}
with open(os.path.join(OUT, "lowk_investigation.json"), "w") as f:
    json.dump(out, f, indent=2)
print(f"\nsaved: results/T2_NG15yr/covariance/lowk_investigation.json")

# ── figure ─────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ax = axes[0]
for r in results:
    n = r["n_bins_used"]
    # recompute eigenspectrum for plotting
    sub = rho[:, :n]
    sigma = np.cov(sub, rowvar=False)
    d = np.sqrt(np.diag(sigma))
    corr = sigma / np.outer(d, d)
    ev = np.sort(np.linalg.eigvalsh(corr))[::-1]
    ax.plot(np.arange(1, n+1), ev, "o-", label=f"n={n}")
ax.axhline(1.0, color="k", lw=0.5, ls="--")
ax.set_xlabel("mode index"); ax.set_ylabel("correlation-matrix eigenvalue")
ax.set_title("(a) Eigenspectrum vs. #bins retained")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

ax = axes[1]
ns = [r["n_bins_used"] for r in results]
ax.plot(ns, [r["eigval_max"] for r in results], "o-", label="λ_max")
ax.plot(ns, [r["max_abs_offdiag_corr"] for r in results], "s-", label="max |off-diag corr|")
ax.plot(ns, [r["mean_abs_offdiag_corr"] for r in results], "^-", label="mean |off-diag corr|")
ax.axhline(0.3, ls=":", color="grey", label="hard threshold")
ax.set_xlabel("#bins retained (low-k first)"); ax.set_ylabel("correlation metric")
ax.set_title("(b) Correlation strength vs. subspace size")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(OUT, "lowk_eigenspectrum.png"), dpi=130)
print(f"saved figure: results/T2_NG15yr/covariance/lowk_eigenspectrum.png")
