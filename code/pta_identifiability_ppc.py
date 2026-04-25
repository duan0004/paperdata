#!/usr/bin/env python3
"""P6 identifiability and posterior-summary predictive diagnostics.

The diagnostics here deliberately stay at the public posterior-summary bridge
level.  They do not claim to be full timing-residual posterior predictive
checks.
"""

from __future__ import annotations

import csv
import math
import sys
from datetime import datetime
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BRIDGE = ROOT / "results" / "prl_reference_bridge"
IDENT = ROOT / "results" / "identifiability"
PPC = ROOT / "results" / "ppc"

if str(ROOT / "code") not in sys.path:
    sys.path.insert(0, str(ROOT / "code"))

from prl_reference_bridge_pipeline import (  # noqa: E402
    DEFAULT_N_FREQ,
    EPTAFreeSpectrumLUT,
    NG15FreeSpectrumLUT,
    PPTAFreeSpectrumLUT,
    bridge_models,
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def fmt(x: float, ndigits: int = 3) -> str:
    if not math.isfinite(float(x)):
        return "nan"
    return f"{float(x):.{ndigits}f}"


def frequency_cut_summary() -> list[dict]:
    rows = [r for r in read_csv(BRIDGE / "frequency_cut_evidence.csv") if r.get("status") == "pass"]
    out: list[dict] = []
    order = ["le_5_nHz", "le_8_nHz", "le_12_nHz", "le_20_nHz", "all"]
    for cut in order:
        cut_rows = [r for r in rows if r["fmax_nHz"] == cut]
        if not cut_rows:
            continue
        best = max(cut_rows, key=lambda r: float(r["ln_Z"]))
        for family in ["SIGW-like", "Cosmo-other", "Astro-curved"]:
            fam_rows = [r for r in cut_rows if r["family"] == family]
            if not fam_rows:
                continue
            fam_best = max(fam_rows, key=lambda r: float(r["ln_Z"]))
            out.append(
                {
                    "fmax_nHz": cut,
                    "best_overall_model": best["model"],
                    "best_overall_family": best["family"],
                    "family": family,
                    "best_family_model": fam_best["model"],
                    "delta_family_best_to_overall": float(fam_best["ln_Z"]) - float(best["ln_Z"]),
                    "n_ppta": best.get("n_ppta", ""),
                    "n_ng15": best.get("n_ng15", ""),
                    "n_epta": best.get("n_epta", ""),
                    "evidence_scale": "hybrid3 anchored posterior-summary bridge",
                }
            )
    return out


def family_gap_summary() -> tuple[list[dict], dict[str, float]]:
    fam_rows = [
        r for r in read_csv(BRIDGE / "family_evidence.csv")
        if r.get("tier") == "hybrid3" and r.get("mode") == "all"
    ]
    if not fam_rows:
        return [], {}
    best = max(float(r["ln_Z_family"]) for r in fam_rows)
    out = []
    gap_by_family: dict[str, float] = {}
    for r in sorted(fam_rows, key=lambda x: float(x["ln_Z_family"]), reverse=True):
        gap = float(r["ln_Z_family"]) - best
        gap_by_family[r["family"]] = gap
        out.append(
            {
                "family": r["family"],
                "n_models": r["n_models"],
                "ln_Z_family": r["ln_Z_family"],
                "delta_to_best_family": gap,
                "interpretation": "equal-weight tested-representative family mixture; not full source-class marginalization",
            }
        )
    return out, gap_by_family


def gamma_projection_summary() -> list[dict]:
    rows = read_csv(BRIDGE / "gamma_projection.csv")
    grouped: dict[str, list[float]] = {}
    losses: dict[str, list[float]] = {}
    for r in rows:
        grouped.setdefault(r["model"], []).append(float(r["gamma_eff_median"]))
        losses.setdefault(r["model"], []).append(float(r["projection_loss_median"]))
    out = []
    for model, vals in sorted(grouped.items()):
        out.append(
            {
                "model": model,
                "gamma_eff_min_median": min(vals),
                "gamma_eff_max_median": max(vals),
                "gamma_eff_span": max(vals) - min(vals),
                "projection_loss_median_max": max(losses[model]),
                "diagnostic_scope": "effective power-law projection across public PTA frequency windows",
            }
        )
    return out


def load_likelihoods():
    return {
        "PPTA-DR3": PPTAFreeSpectrumLUT(n_freq=DEFAULT_N_FREQ),
        "NANOGrav-15yr": NG15FreeSpectrumLUT(n_freq=DEFAULT_N_FREQ, analysis="hd"),
        "EPTA-DR2new": EPTAFreeSpectrumLUT(n_freq=DEFAULT_N_FREQ, dataset="DR2new"),
    }


def model_ppc_rows() -> tuple[list[dict], list[dict]]:
    models = bridge_models()
    lnz_rows = [r for r in read_csv(BRIDGE / "hybrid3_bridge_lnZ.csv") if r.get("status") == "pass"]
    lnz_rows.sort(key=lambda r: float(r["ln_Z"]), reverse=True)
    selected = [r["model"] for r in lnz_rows[:7]]
    likes = load_likelihoods()
    per_bin: list[dict] = []
    summary: list[dict] = []

    data_summaries = {name: like.summary_per_bin() for name, like in likes.items()}
    for key in selected:
        sample_path = BRIDGE / "P2_hybrid3" / key / "posterior_equal_weight.npz"
        if not sample_path.exists():
            continue
        model = models[key]
        data = np.load(sample_path, allow_pickle=True)
        samples = np.asarray(data["samples"], dtype=float)
        pnames = [str(x) for x in data["param_names"]]
        idx = np.linspace(0, samples.shape[0] - 1, min(512, samples.shape[0]), dtype=int)
        for dname, like in likes.items():
            pred = []
            for ii in idx:
                params = dict(model.param_defaults)
                for j, pname in enumerate(pnames):
                    params[pname] = float(samples[ii, j])
                pred.append(like.model_to_log10_rho(model, params))
            pred_arr = np.asarray(pred)
            low_abs_z = []
            low_covered = []
            for k in range(like.n_freq):
                p16, p50, p84 = np.percentile(pred_arr[:, k], [16, 50, 84])
                d = data_summaries[dname][k]
                sigma = max(float(d["log10rho_std"]), 1e-6)
                z = (float(p50) - float(d["log10rho_mean"])) / sigma
                covered = float(p16) <= float(d["log10rho_mean"]) <= float(p84)
                if k < 4:
                    low_abs_z.append(abs(z))
                    low_covered.append(covered)
                per_bin.append(
                    {
                        "model": key,
                        "family": getattr(model, "family", ""),
                        "dataset": dname,
                        "bin": k + 1,
                        "f_nHz": float(like.f_bins[k] * 1e9),
                        "data_log10rho_mean": float(d["log10rho_mean"]),
                        "data_log10rho_std": float(d["log10rho_std"]),
                        "pred_log10rho_p16": float(p16),
                        "pred_log10rho_median": float(p50),
                        "pred_log10rho_p84": float(p84),
                        "z_median": float(z),
                        "data_mean_inside_model_68": int(covered),
                        "diagnostic_scope": "posterior-summary residual diagnostic, not full timing-residual PPC",
                    }
                )
            summary.append(
                {
                    "model": key,
                    "family": getattr(model, "family", ""),
                    "dataset": dname,
                    "n_low_bins": 4,
                    "mean_abs_z_low4": float(np.mean(low_abs_z)),
                    "max_abs_z_low4": float(np.max(low_abs_z)),
                    "coverage_low4": float(np.mean(low_covered)),
                    "n_projected_samples": int(len(idx)),
                }
            )
    return per_bin, summary


def write_report(
    cut_rows: list[dict],
    family_rows: list[dict],
    gamma_rows: list[dict],
    ppc_summary: list[dict],
) -> None:
    best_family_gap = {r["family"]: float(r["delta_to_best_family"]) for r in family_rows}
    astro_gap = best_family_gap.get("Astro-curved", float("nan"))
    gamma_min = min(float(r["gamma_eff_min_median"]) for r in gamma_rows) if gamma_rows else float("nan")
    gamma_max = max(float(r["gamma_eff_max_median"]) for r in gamma_rows) if gamma_rows else float("nan")

    cut_all = [r for r in cut_rows if r["fmax_nHz"] == "all" and r["family"] == "Astro-curved"]
    astro_all_gap = float(cut_all[0]["delta_family_best_to_overall"]) if cut_all else float("nan")

    low4 = [r for r in ppc_summary if r["family"] in {"Astro-curved", "SIGW-like", "Cosmo-other"}]
    mean_low_z = float(np.mean([float(r["mean_abs_z_low4"]) for r in low4])) if low4 else float("nan")

    lines = [
        "# P6 Identifiability and Posterior-Summary Predictive Diagnostics",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Scope",
        "",
        "These diagnostics use the hybrid3 anchored posterior-summary bridge and public per-bin free-spectrum likelihood summaries. They are not a full timing-residual posterior predictive check and must not be cited as a full IPTA-level likelihood result.",
        "",
        "## Identifiability",
        "",
        f"- Equal-weight tested-representative family mixture: Astro-curved is {fmt(astro_gap)} nat below the leading family on the hybrid3 bridge scale.",
        f"- In the full 14-bin frequency-cut diagnostic, the best Astro-curved template is {fmt(astro_all_gap)} nat below the best tested template.",
        f"- Curved-SMBHB effective-slope projections span gamma_eff median values from {fmt(gamma_min)} to {fmt(gamma_max)} across the three public PTA frequency windows.",
        "",
        "Interpretation: the present bridge diagnostics identify an amplitude-like mode and a low-frequency-curvature mode more clearly than a unique source class. Higher-order source-family structure remains order-unity in log evidence for the tested representatives.",
        "",
        "## Low-Frequency Posterior-Summary Residuals",
        "",
        f"- Across selected leading templates and the three public posterior-summary likelihoods, the mean low-four-bin absolute residual is {fmt(mean_low_z)} in LUT-standard-deviation units.",
        "- Residual tables are written per model, dataset, and bin so the statement can be checked without relying on a plotted summary.",
        "",
        "## Files",
        "",
        f"- `{(IDENT / 'low_frequency_ranking.csv').relative_to(ROOT)}`",
        f"- `{(IDENT / 'family_gap_summary.csv').relative_to(ROOT)}`",
        f"- `{(IDENT / 'gamma_projection_summary.csv').relative_to(ROOT)}`",
        f"- `{(IDENT / 'spectral_basis_evidence.csv').relative_to(ROOT)}`",
        f"- `{(PPC / 'posterior_summary_residuals.csv').relative_to(ROOT)}`",
        f"- `{(PPC / 'family_ppc_lowfreq.csv').relative_to(ROOT)}`",
    ]
    (IDENT / "mode_compression.md").write_text("\n".join(lines) + "\n")


def write_basis_evidence(family_rows: list[dict], gamma_rows: list[dict]) -> None:
    family_best = {r["family"]: float(r["delta_to_best_family"]) for r in family_rows}
    gamma_span = (
        max(float(r["gamma_eff_max_median"]) for r in gamma_rows) -
        min(float(r["gamma_eff_min_median"]) for r in gamma_rows)
    ) if gamma_rows else float("nan")
    rows = [
        {
            "mode": "amplitude_like_common_power",
            "diagnostic": "finite timing-likelihood pilot plus bridge posterior samples",
            "value": "finite pilot likelihood; bridge posteriors available",
            "interpretation": "baseline amplitude mode is technically measurable but not yet production timing evidence",
        },
        {
            "mode": "low_frequency_curvature",
            "diagnostic": "gamma_eff projection across PPTA/NG15/EPTA windows",
            "value": f"gamma_eff median span {gamma_span:.6g}",
            "interpretation": "curved templates project to steep effective slopes in the public PTA windows",
        },
        {
            "mode": "source_family_label",
            "diagnostic": "hybrid3 tested-representative family mixture",
            "value": f"Astro-curved delta {family_best.get('Astro-curved', float('nan')):.6g} nat",
            "interpretation": "source-family separation remains order-unity for the tested representatives",
        },
    ]
    write_csv(
        IDENT / "spectral_basis_evidence.csv",
        rows,
        ["mode", "diagnostic", "value", "interpretation"],
    )


def write_figure(per_bin: list[dict]) -> None:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except Exception as exc:
        (PPC / "family_ppc_figures.skipped.txt").write_text(f"matplotlib unavailable: {exc}\n")
        return

    PPC.mkdir(parents=True, exist_ok=True)
    datasets = ["PPTA-DR3", "NANOGrav-15yr", "EPTA-DR2new"]
    models = []
    for row in per_bin:
        if row["model"] not in models:
            models.append(row["model"])
    with PdfPages(PPC / "family_ppc_figures.pdf") as pdf:
        for dname in datasets:
            fig, ax = plt.subplots(figsize=(8, 4.5))
            for model in models[:7]:
                rows = [r for r in per_bin if r["dataset"] == dname and r["model"] == model and int(r["bin"]) <= 4]
                if not rows:
                    continue
                x = [int(r["bin"]) for r in rows]
                y = [float(r["z_median"]) for r in rows]
                ax.plot(x, y, marker="o", linewidth=1.2, label=model)
            ax.axhline(0.0, color="black", linewidth=0.8)
            ax.axhspan(-1.0, 1.0, color="0.9", zorder=-1)
            ax.set_title(f"Low-frequency posterior-summary residuals: {dname}")
            ax.set_xlabel("frequency bin")
            ax.set_ylabel("median residual / LUT std")
            ax.set_xticks([1, 2, 3, 4])
            ax.legend(fontsize=7, ncol=2)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def main() -> int:
    IDENT.mkdir(parents=True, exist_ok=True)
    PPC.mkdir(parents=True, exist_ok=True)

    cut_rows = frequency_cut_summary()
    family_rows, _ = family_gap_summary()
    gamma_rows = gamma_projection_summary()
    per_bin, ppc_summary = model_ppc_rows()

    write_csv(
        IDENT / "low_frequency_ranking.csv",
        cut_rows,
        [
            "fmax_nHz",
            "best_overall_model",
            "best_overall_family",
            "family",
            "best_family_model",
            "delta_family_best_to_overall",
            "n_ppta",
            "n_ng15",
            "n_epta",
            "evidence_scale",
        ],
    )
    write_csv(
        IDENT / "family_gap_summary.csv",
        family_rows,
        ["family", "n_models", "ln_Z_family", "delta_to_best_family", "interpretation"],
    )
    write_csv(
        IDENT / "gamma_projection_summary.csv",
        gamma_rows,
        [
            "model",
            "gamma_eff_min_median",
            "gamma_eff_max_median",
            "gamma_eff_span",
            "projection_loss_median_max",
            "diagnostic_scope",
        ],
    )
    write_basis_evidence(family_rows, gamma_rows)
    write_csv(
        PPC / "posterior_summary_residuals.csv",
        per_bin,
        [
            "model",
            "family",
            "dataset",
            "bin",
            "f_nHz",
            "data_log10rho_mean",
            "data_log10rho_std",
            "pred_log10rho_p16",
            "pred_log10rho_median",
            "pred_log10rho_p84",
            "z_median",
            "data_mean_inside_model_68",
            "diagnostic_scope",
        ],
    )
    write_csv(
        PPC / "family_ppc_lowfreq.csv",
        ppc_summary,
        [
            "model",
            "family",
            "dataset",
            "n_low_bins",
            "mean_abs_z_low4",
            "max_abs_z_low4",
            "coverage_low4",
            "n_projected_samples",
        ],
    )
    write_figure(per_bin)
    write_report(cut_rows, family_rows, gamma_rows, ppc_summary)
    print({"cut_rows": len(cut_rows), "ppc_bins": len(per_bin), "ppc_summary": len(ppc_summary)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
