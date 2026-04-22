#!/usr/bin/env python3
"""
template_projection.py — Plan D / T1b.3
=========================================

Project the three Bayes-factor-refit templates (SMBHB, SIGW-Gauss,
cs_stable) at their posterior-mean parameters onto the eigenmodes of
the 14-bin posterior covariance Σ_bin computed in T1b.1.

Discriminating metric D_M:
    D_M = Σ_k λ_k^{-1} (c_k^M - c_k^SMBHB)^2
where c_k^M = e_k · (μ_template_M - μ_data) / stdev_bin.
λ_k^{-1} weights mean "low-variance direction = strong discriminator".

Expected magnitude: given max eigenvalue ≈ 1.48 and min ≈ 0.71 in the
correlation matrix, the D metric per bin-pair is O(1).  Total discriminating
power is ~1 nat for each 1 σ separation between models in the dominant
eigendirection.  Cf. pilot memo 2026-04-21: this is a *modest* observable.

Outputs:
  results/T2_NG15yr/covariance/template_projection.json
  results/T2_NG15yr/covariance/template_projection.png
"""
import os, sys, json, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import gwb_templates as gwt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COV  = os.path.join(BASE, "results", "T2_NG15yr", "covariance")
BF   = os.path.join(BASE, "results", "T2_NG15yr", "bayes_factors", "bayes_factors.json")

# Load eigenstructure
mu_bin     = np.load(os.path.join(COV, "mu_bin.npy"))
sigma_bin  = np.load(os.path.join(COV, "sigma_bin.npy"))
corr_bin   = np.load(os.path.join(COV, "corr_bin.npy"))
eigvals    = np.load(os.path.join(COV, "eigvals.npy"))
eigvecs    = np.load(os.path.join(COV, "eigvecs.npy"))
d_bin      = np.sqrt(np.diag(sigma_bin))

# Load Bayes-factor fits
with open(BF) as f:
    bf = json.load(f)
FREQS    = np.array(bf["freqs_Hz"])
TSPAN_S  = bf["Tspan_yr"] * 365.25 * 86400
H0       = gwt.H0
N_BINS   = len(FREQS)

def omega_to_log10_rho(omega, f):
    omega = np.asarray(omega)
    omega = np.clip(omega, 1e-40, np.inf)
    rho2 = H0**2 * omega / (8.0 * np.pi**4 * f**5 * TSPAN_S)
    return 0.5 * np.log10(rho2)

# ── posterior-mean template spectra ─────────────────────────────────
smbhb = gwt.get_template('smbhb')
sigw  = gwt.get_template('sigw_gauss')
cs    = gwt.get_template('cosmic_strings_stable')

# SMBHB: single param [log10_A], gamma fixed to 13/3
log10_A_smbhb = bf["results"]["smbhb"]["posterior_mean"][0]
spec_smbhb = omega_to_log10_rho(
    smbhb.omega_gw(FREQS, log10_A=log10_A_smbhb, gamma=13.0/3.0),
    FREQS,
)

# SIGW-Gauss: [log10_As, log10_fstar, Delta]
pm_sigw = bf["results"]["sigw_gauss"]["posterior_mean"]
spec_sigw = omega_to_log10_rho(
    sigw.omega_gw(FREQS, log10_As=pm_sigw[0],
                  log10_f_star=pm_sigw[1], Delta=pm_sigw[2]),
    FREQS,
)

# CS-stable: [log10_Gmu]
log10_Gmu_cs = bf["results"]["cs_stable"]["posterior_mean"][0]
spec_cs = omega_to_log10_rho(
    cs.omega_gw(FREQS, log10_Gmu=log10_Gmu_cs), FREQS,
)

print(f"SMBHB       (log10_A = {log10_A_smbhb:+.3f}):")
print(f"  pred log10_ρ = {np.array2string(spec_smbhb, precision=2)}")
print(f"  data  μ_bin  = {np.array2string(mu_bin, precision=2)}")
print(f"  residual     = {np.array2string(spec_smbhb - mu_bin, precision=2)}")
print()
print(f"SIGW-Gauss  (log10_As={pm_sigw[0]:+.3f}, log10_f*={pm_sigw[1]:+.3f}, Δ={pm_sigw[2]:.3f}):")
print(f"  pred log10_ρ = {np.array2string(spec_sigw, precision=2)}")
print()
print(f"CS-stable   (log10_Gμ = {log10_Gmu_cs:+.3f}):")
print(f"  pred log10_ρ = {np.array2string(spec_cs, precision=2)}")

# ── whiten residuals + project ──────────────────────────────────────
# whitened residual r_M = (spec_M - μ_bin) / d_bin   (14-vector)
# projection coeffs c_k^M = e_k · r_M   (where e_k are columns of eigvecs)
def project(spec):
    r = (spec - mu_bin) / d_bin
    c = eigvecs.T @ r    # shape (14,)
    return r, c

r_smbhb, c_smbhb = project(spec_smbhb)
r_sigw,  c_sigw  = project(spec_sigw)
r_cs,    c_cs    = project(spec_cs)

print("\n── Projections onto correlation-matrix eigenmodes ──")
hdr = "mode    λ_k    c_SMBHB    c_SIGW     c_CS     "
print(hdr); print("-" * len(hdr))
for k in range(N_BINS):
    print(f"{k+1:2d}   {eigvals[k]:.3f}  {c_smbhb[k]:+7.3f}   {c_sigw[k]:+7.3f}   {c_cs[k]:+7.3f}")

# ── Mahalanobis-style D metric vs SMBHB baseline ────────────────────
# D_M = Σ_k λ_k^{-1} (c_k^M - c_k^SMBHB)^2
def D_metric(c_M):
    diff = c_M - c_smbhb
    return float(np.sum(diff**2 / eigvals))

D_sigw = D_metric(c_sigw)
D_cs   = D_metric(c_cs)

# Also compute Δ² = (r_M - r_smbhb)^T Σ^{-1} (r_M - r_smbhb) for reference
# where Σ here is the correlation matrix (eigenbasis). This is CAR minus ceffyl:
#   CAR uses full Σ^{-1} = Σ_k λ^{-1} e_k e_k^T
#   ceffyl approximates Σ^{-1} ≈ I (diagonal in original basis)
# so the Δ between the two metrics quantifies the bin-correlation effect.
def chi2_diag(spec_M):
    # ceffyl-style: diag-only, weight 1/1 = 1
    return float(np.sum((spec_M - spec_smbhb)**2 / d_bin**2))

def chi2_CAR(spec_M):
    # CAR-style: full correlation inverse
    r = (spec_M - spec_smbhb) / d_bin
    Cinv = eigvecs @ np.diag(1.0/eigvals) @ eigvecs.T
    return float(r @ Cinv @ r)

chi2_sigw_diag = chi2_diag(spec_sigw)
chi2_sigw_CAR  = chi2_CAR(spec_sigw)
chi2_cs_diag   = chi2_diag(spec_cs)
chi2_cs_CAR    = chi2_CAR(spec_cs)

print(f"\n── Discriminating metrics (vs SMBHB at its posterior mean) ──")
print(f"  SIGW-Gauss:  D = {D_sigw:.3f}   χ²_diag = {chi2_sigw_diag:.3f}   χ²_CAR = {chi2_sigw_CAR:.3f}"
      f"   Δχ² = {chi2_sigw_CAR - chi2_sigw_diag:+.3f}")
print(f"  CS-stable :  D = {D_cs:.3f}   χ²_diag = {chi2_cs_diag:.3f}   χ²_CAR = {chi2_cs_CAR:.3f}"
      f"   Δχ² = {chi2_cs_CAR - chi2_cs_diag:+.3f}")

print(f"\n  Δχ² is the magnitude of the CAR vs ceffyl correction per template,")
print(f"  and is a direct proxy for the ΔlnZ = -½ Δχ² that CAR would introduce.")
print(f"  Half-Δχ² for SIGW:  {-0.5*(chi2_sigw_CAR-chi2_sigw_diag):+.3f} nats")
print(f"  Half-Δχ² for CS:    {-0.5*(chi2_cs_CAR-chi2_cs_diag):+.3f} nats")

# ── save ────────────────────────────────────────────────────────────
out = {
    "task": "T1b.3",
    "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
    "templates_evaluated_at": {
        "smbhb":      {"log10_A": log10_A_smbhb, "gamma": 13.0/3.0},
        "sigw_gauss": {"log10_As": pm_sigw[0], "log10_fstar": pm_sigw[1], "Delta": pm_sigw[2]},
        "cs_stable":  {"log10_Gmu": log10_Gmu_cs},
    },
    "spec_log10_rho": {
        "smbhb":      spec_smbhb.tolist(),
        "sigw_gauss": spec_sigw.tolist(),
        "cs_stable":  spec_cs.tolist(),
    },
    "projections_onto_corr_eigenmodes": {
        "smbhb":      c_smbhb.tolist(),
        "sigw_gauss": c_sigw.tolist(),
        "cs_stable":  c_cs.tolist(),
    },
    "D_metric_vs_smbhb": {
        "sigw_gauss": D_sigw,
        "cs_stable":  D_cs,
    },
    "chi2_comparison": {
        "sigw_gauss": {"diag": chi2_sigw_diag, "CAR": chi2_sigw_CAR,
                       "delta_chi2_CAR_minus_diag": chi2_sigw_CAR - chi2_sigw_diag,
                       "implied_delta_lnZ_CAR_correction": -0.5*(chi2_sigw_CAR - chi2_sigw_diag)},
        "cs_stable":  {"diag": chi2_cs_diag, "CAR": chi2_cs_CAR,
                       "delta_chi2_CAR_minus_diag": chi2_cs_CAR - chi2_cs_diag,
                       "implied_delta_lnZ_CAR_correction": -0.5*(chi2_cs_CAR - chi2_cs_diag)},
    },
}
with open(os.path.join(COV, "template_projection.json"), "w") as f:
    json.dump(out, f, indent=2)

# ── figure ─────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

ax = axes[0]
idx = np.arange(1, N_BINS+1)
ax.errorbar(idx, mu_bin, yerr=d_bin, fmt="o", color="k",
            label="data μ ± σ", capsize=2, ms=4)
ax.plot(idx, spec_smbhb, "s-", color="C0", label="SMBHB (post mean)")
ax.plot(idx, spec_sigw,  "^-", color="C1", label="SIGW-Gauss (post mean)")
ax.plot(idx, spec_cs,    "v-", color="C2", label="CS-stable (post mean)")
ax.set_xlabel("bin index"); ax.set_ylabel("log10 ρ")
ax.set_title("(a) Template predictions vs NG15yr HD free-spectrum posterior")
ax.legend(fontsize=9); ax.grid(alpha=0.3)

ax = axes[1]
ax.plot(idx, c_smbhb, "s-", color="C0", label=f"SMBHB")
ax.plot(idx, c_sigw,  "^-", color="C1", label=f"SIGW-Gauss (D={D_sigw:.2f})")
ax.plot(idx, c_cs,    "v-", color="C2", label=f"CS-stable (D={D_cs:.2f})")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("correlation-matrix mode index k")
ax.set_ylabel("projection coeff $c_k$")
ax.set_title("(b) Template projections onto Σ_bin eigenmodes")
ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(COV, "template_projection.png"), dpi=130)
print(f"\nsaved: {os.path.relpath(os.path.join(COV, 'template_projection.json'), BASE)}")
print(f"saved: {os.path.relpath(os.path.join(COV, 'template_projection.png'), BASE)}")
