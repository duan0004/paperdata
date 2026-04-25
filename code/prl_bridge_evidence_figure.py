#!/usr/bin/env python3
"""Build the PRL cross-PTA bridge evidence figure.

The figure is a visualization layer only.  It reads the production P0--P6
outputs written by code/prl_reference_bridge_pipeline.py and does not alter
the evidence calculations.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
IN_DIR = ROOT / "results" / "prl_reference_bridge"
OUT_DIR = ROOT / "results" / "T2_NG15yr" / "figures"

LOCAL3 = IN_DIR / "local3_bridge_lnZ.csv"
HYBRID3 = IN_DIR / "hybrid3_bridge_lnZ.csv"
FAMILY = IN_DIR / "family_evidence.csv"

OUT_PDF = OUT_DIR / "prl_bridge_evidence_figure.pdf"
OUT_PNG = OUT_DIR / "prl_bridge_evidence_figure.png"

FAMILY_COLORS = {
    "SIGW-like": "#31688e",
    "Cosmo-other": "#35b779",
    "Astro-curved": "#f89540",
    "Astro-simple": "#7a7a7a",
}

DISPLAY = {
    "Cosmic-Superstrings": "Cosmic superstrings",
    "SIGW-Delta": "SIGW-delta",
    "PBH-SIGW-Analytic": "PBH-SIGW analytic",
    "SIGW-Gaussian": "SIGW-Gaussian",
    "SMBHB-Eccentric": "Ecc.-curved SMBHB",
    "SMBHB-BrokenPL": "Broken-PL SMBHB",
    "SMBHB-Env": "Env. SMBHB",
    "SMBHB-Turnover": "Turnover SMBHB",
    "SMBHB-PowerLaw": "Power-law SMBHB",
}


def _read_rank(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df["status"].eq("pass")].copy()
    best_err = float(df.loc[df["ln_Z"].idxmax(), "ln_Z_err"])
    df["delta_pos"] = df["delta_to_best"].astype(float)
    df["delta_err"] = (df["ln_Z_err"].astype(float) ** 2 + best_err**2) ** 0.5
    df["label"] = df["model"].map(DISPLAY).fillna(df["model"])
    return df.sort_values("delta_to_best", ascending=False)


def _family_table() -> pd.DataFrame:
    df = pd.read_csv(FAMILY)
    df = df[(df["tier"].eq("hybrid3")) & (df["mode"].eq("all"))].copy()
    best = df["ln_Z_family"].max()
    df["delta"] = df["ln_Z_family"] - best
    order = ["SIGW-like", "Astro-curved", "Cosmo-other", "Astro-simple"]
    df["family"] = pd.Categorical(df["family"], order, ordered=True)
    return df.sort_values("family")


def build() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    local = _read_rank(LOCAL3)
    hybrid = _read_rank(HYBRID3)
    top_models = [
        "Cosmic-Superstrings",
        "SIGW-Delta",
        "PBH-SIGW-Analytic",
        "SIGW-Gaussian",
        "SMBHB-Eccentric",
        "SMBHB-BrokenPL",
        "SMBHB-Env",
    ]
    local = local.set_index("model").loc[top_models].reset_index()
    hybrid = hybrid.set_index("model").loc[top_models].reset_index()
    families = _family_table()

    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.35), gridspec_kw={"width_ratios": [1.35, 1.0]})
    ax = axes[0]
    y = list(range(len(top_models)))
    colors = [FAMILY_COLORS[f] for f in hybrid["family"]]

    for yi, row, color in zip(y, local.itertuples(), colors):
        ax.errorbar(
            row.delta_pos,
            yi,
            xerr=row.delta_err,
            fmt="o",
            ms=4.8,
            mfc="white",
            mec=color,
            mew=1.2,
            ecolor=color,
            elinewidth=0.65,
            capsize=1.8,
            label="local3" if yi == y[0] else None,
        )
    for yi, row, color in zip(y, hybrid.itertuples(), colors):
        ax.errorbar(
            row.delta_pos,
            yi,
            xerr=row.delta_err,
            fmt="s",
            ms=4.8,
            mfc=color,
            mec="black",
            mew=0.35,
            ecolor=color,
            elinewidth=0.65,
            capsize=1.8,
            label="hybrid3" if yi == y[0] else None,
        )
    for yi, xl, xh in zip(y, local["delta_pos"], hybrid["delta_pos"]):
        ax.plot([xl, xh], [yi, yi], color="#b8b8b8", lw=0.8, zorder=0)

    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(-1, color="#666666", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(hybrid["label"])
    ax.invert_yaxis()
    ax.set_xlim(-1.7, 0.12)
    ax.set_xlabel(r"$\Delta\ln Z$ relative to best model")
    ax.set_title("A. Cross-PTA bridge ranking")
    ax.grid(axis="x", color="#dedede", linewidth=0.5)
    ax.legend(loc="lower left", frameon=False, handletextpad=0.35)

    ax2 = axes[1]
    fam_colors = [FAMILY_COLORS[str(f)] for f in families["family"]]
    ax2.barh(range(len(families)), families["delta"], color=fam_colors, edgecolor="black", linewidth=0.35)
    ax2.set_yticks(range(len(families)))
    ax2.set_yticklabels(families["family"].astype(str))
    ax2.invert_yaxis()
    ax2.axvline(0, color="black", lw=0.8)
    ax2.axvline(-1, color="#666666", lw=0.8, ls="--")
    ax2.set_xlim(-1.25, 0.08)
    ax2.set_xlabel(r"$\Delta\ln Z_\mathrm{family}$")
    ax2.set_title("B. Hybrid3 family evidence")
    ax2.grid(axis="x", color="#dedede", linewidth=0.5)
    for idx, val in enumerate(families["delta"]):
        ax2.text(-0.08, idx, f"{val:.2f}", ha="right", va="center", color="black")

    fig.tight_layout(pad=1.05)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=220, bbox_inches="tight")
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    build()
