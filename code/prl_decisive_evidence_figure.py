#!/usr/bin/env python3
"""Build the PRL decisive evidence figure.

Panel A shows the evidence-ablation path for SIGW-Gaussian.  Panel B compares
cosmological stochastic templates with the curved-SMBHB control family
on the matched official PTArcade ceffyl-density evidence scale.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BFDIR = ROOT / "results/T2_NG15yr/bayes_factors"
FIGDIR = ROOT / "results/T2_NG15yr/figures"
OUT_PDF = FIGDIR / "prl_decisive_evidence_figure.pdf"
OUT_PNG = FIGDIR / "prl_decisive_evidence_figure.png"
OUT_JSON = FIGDIR / "prl_decisive_evidence_figure_data.json"
OUT_MD = FIGDIR / "prl_decisive_evidence_figure.md"

ROBUST_JSON = BFDIR / "prl_official_density_robustness.json"
H7_JSON = BFDIR / "prl_H7_astrophysical_family.json"
ABLATION_MD = BFDIR / "prl_ablation_table.md"


def load_json(path: Path):
    return json.loads(path.read_text())


def combined_uncertainty(sample_std, nested_err):
    return float(np.sqrt(float(sample_std) ** 2 + float(nested_err) ** 2))


def build_data():
    robust = load_json(ROBUST_JSON)
    h7 = load_json(H7_JSON)

    # Values are from results/T2_NG15yr/bayes_factors/prl_ablation_table.md.
    ablation = [
        {
            "label": "local template\nlocal KDE",
            "short": "local-template/local-KDE",
            "baseline": "fixed-gamma SMBHB",
            "lnB": -1.590,
            "err": 0.154,
            "residual_vs_published": -5.633,
        },
        {
            "label": "official template\nlocal KDE",
            "short": "official-template/local-KDE",
            "baseline": "PTArcade BHB prior",
            "lnB": 0.687,
            "err": 0.147,
            "residual_vs_published": -3.356,
        },
        {
            "label": "official template\nofficial density\nseed 42",
            "short": "official-template/official-density seed42",
            "baseline": "PTArcade BHB prior",
            "lnB": 4.225,
            "err": 0.173,
            "residual_vs_published": 0.182,
        },
        {
            "label": "official template\nofficial density\n5-config mean",
            "short": "official-template/official-density robust",
            "baseline": "PTArcade BHB prior",
            "lnB": 4.520,
            "err": 0.185,
            "residual_vs_published": 0.477,
        },
    ]

    cosmology_labels = {
        "sigw_gauss_ptarcade": "SIGW Gaussian",
        "sigw_delta_ptarcade": "SIGW delta",
        "super_ptarcade": "Cosmic superstrings",
    }
    model_rows = []
    for key, label in cosmology_labels.items():
        item = robust["aggregate"][key]
        model_rows.append(
            {
                "group": "cosmology",
                "model": key,
                "label": label,
                "lnB": item["lnB_mean"],
                "sample_std": item["lnB_sample_std"],
                "nested_err": item["mean_nested_lnB_err"],
                "err": combined_uncertainty(
                    item["lnB_sample_std"], item["mean_nested_lnB_err"]
                ),
            }
        )

    astro_labels = {
        "ecc_supp_fixed_gamma": "Eccentricity-inspired\ncurvature",
        "broken_pl_fixed_gamma": "Broken power law\ncurvature",
        "env_fixed_gamma": "Environmental\nturnover",
        "env_free_gamma": "Environmental\nfree gamma",
    }
    for key, label in astro_labels.items():
        item = h7["aggregate"][key]
        model_rows.append(
            {
                "group": "curved SMBHB control",
                "class": item["class"],
                "model": key,
                "label": label,
                "lnB": item["lnB_mean"],
                "sample_std": item["lnB_sample_std"],
                "nested_err": item["mean_nested_lnB_err"],
                "err": combined_uncertainty(
                    item["lnB_sample_std"], item["mean_nested_lnB_err"]
                ),
                "qmc_ti_lnB": (
                    h7["qmc_ti_crosscheck"]["aggregate"].get(key, {}).get("lnB_ti_mean")
                ),
            }
        )

    model_rows.sort(key=lambda row: row["lnB"], reverse=True)
    return {
        "source_files": {
            "ablation": str(ABLATION_MD.relative_to(ROOT)),
            "official_density_robustness": str(ROBUST_JSON.relative_to(ROOT)),
            "astrophysical_family": str(H7_JSON.relative_to(ROOT)),
        },
        "ablation": ablation,
        "model_rows": model_rows,
        "notes": {
            "baseline": "All Panel B Bayes factors are relative to the matched PTArcade BHB-prior SMBHB baseline.",
            "panel_a_baselines": "Panel A compares Bayes factors relative to the baseline indicated for each row: first bar fixed-gamma SMBHB, all remaining bars PTArcade BHB-prior SMBHB.",
            "uncertainty": "Panel B error bars use quadrature of five-configuration sample scatter and mean nested-sampling lnB error.",
            "h7_family_pass": h7["family_pass"],
            "h7_family_pass_classes": h7["family_pass_classes"],
        },
    }


def make_figure(data):
    plt.rcParams.update(
        {
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.0,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.7,
            "legend.fontsize": 7.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), constrained_layout=True)

    ax = axes[0]
    ablation = data["ablation"]
    x = np.arange(len(ablation))
    vals = np.array([row["lnB"] for row in ablation])
    errs = np.array([row["err"] for row in ablation])
    colors = ["#8d99ae", "#4f6d7a", "#2a9d8f", "#1f7a5c"]
    ax.bar(x, vals, yerr=errs, color=colors, edgecolor="black", linewidth=0.4, capsize=2.5)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axhline(4.04305126783455, color="#7f1d1d", linestyle="--", linewidth=0.9)
    ax.text(
        0.03,
        4.12,
        "published SIGW-Gaussian lnB",
        color="#7f1d1d",
        transform=ax.get_yaxis_transform(),
        va="bottom",
    )
    ax.set_ylabel(r"$\ln B$ vs row-specific baseline")
    ax.set_title("A. Evidence calibration")
    ax.set_xticks(x)
    ax.set_xticklabels([row["label"] for row in ablation], rotation=30, ha="right")
    ax.set_ylim(-2.2, 5.2)
    ax.grid(axis="y", alpha=0.22, linewidth=0.6)

    ax = axes[1]
    rows = data["model_rows"]
    y = np.arange(len(rows))[::-1]
    vals = np.array([row["lnB"] for row in rows])
    errs = np.array([row["err"] for row in rows])
    colors = [
        "#26547c" if row["group"] == "cosmology" else "#b36b00" for row in rows
    ]
    ax.barh(y, vals, xerr=errs, color=colors, edgecolor="black", linewidth=0.4, capsize=2.4)
    ax.axvspan(3.5575775855519254 - 1.0, 4.519926101285057 + 1.0, color="#d9ead3", alpha=0.35)
    ax.axvline(3.5575775855519254, color="#52796f", linestyle=":", linewidth=0.9)
    ax.axvline(4.519926101285057, color="#52796f", linestyle=":", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([row["label"] for row in rows])
    ax.set_xlabel(r"$\ln B$ vs matched BHB prior")
    ax.set_title("B. Matched-scale source ranking")
    ax.set_xlim(2.2, 5.0)
    ax.grid(axis="x", alpha=0.22, linewidth=0.6)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=250)


def write_md(data):
    lines = [
        "# PRL Decisive Evidence Figure",
        "",
        f"**PDF**: `{OUT_PDF.relative_to(ROOT)}`",
        f"**PNG**: `{OUT_PNG.relative_to(ROOT)}`",
        f"**Data JSON**: `{OUT_JSON.relative_to(ROOT)}`",
        "",
        "## Panel A: Calibration Ablation",
        "",
        "| Step | baseline | lnB | err | residual vs published |",
        "|---|---|---:|---:|---:|",
    ]
    for row in data["ablation"]:
        lines.append(
            f"| {row['short']} | {row['baseline']} | `{row['lnB']:+.3f}` | `{row['err']:.3f}` | "
            f"`{row['residual_vs_published']:+.3f}` |"
        )
    lines.extend(
        [
            "",
            "## Panel B: Matched-Scale Ranking",
            "",
            "| Group | Model | lnB | plotted err |",
            "|---|---|---:|---:|",
        ]
    )
    for row in data["model_rows"]:
        lines.append(
            f"| {row['group']} | {row['model']} | `{row['lnB']:+.3f}` | `{row['err']:.3f}` |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The figure shows that the absolute-evidence recovery requires the matched",
            "official template, official density, and matched BHB-prior baseline.  Panel A",
            "uses the row-specific baseline listed above; Panel B uses only the matched",
            "PTArcade BHB-prior baseline.  On that matched scale, multiple curved-SMBHB",
            "controls occupy the same evidence tier as the leading stochastic new-physics",
            "templates.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    data = build_data()
    OUT_JSON.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    make_figure(data)
    write_md(data)
    print(f"saved: {OUT_PDF}")
    print(f"saved: {OUT_PNG}")
    print(f"saved: {OUT_JSON}")
    print(f"saved: {OUT_MD}")


if __name__ == "__main__":
    main()
