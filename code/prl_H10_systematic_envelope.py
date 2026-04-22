#!/usr/bin/env python3
"""PRL H10: assemble the systematic evidence envelope.

The script does not generate new Bayes factors.  It collects completed H1/H6/H7
and H9 outputs into a single systematic budget for the PRL source-identification
claim.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BFDIR = ROOT / "results/T2_NG15yr/bayes_factors"
OUTFILE = BFDIR / "prl_H10_systematic_envelope.json"
OUT_MD = BFDIR / "prl_H10_systematic_envelope.md"

ROBUST_JSON = BFDIR / "prl_official_density_robustness.json"
H456_JSON = BFDIR / "prl_H4_H5_H6_experiments.json"
QMC_JSON = BFDIR / "prl_evidence_ti_qmc_crosscheck.json"
H7_JSON = BFDIR / "prl_H7_astrophysical_family.json"
H9_JSON = BFDIR / "prl_H9_bin_driver_analysis.json"


def load(path: Path):
    return json.loads(path.read_text())


def max_abs(values):
    return float(max(abs(v) for v in values))


def build_envelope():
    robust = load(ROBUST_JSON)
    h456 = load(H456_JSON)
    qmc = load(QMC_JSON)
    h7 = load(H7_JSON)
    h9 = load(H9_JSON)

    cosmology_map = {
        "SIGW-Gaussian": "sigw_gauss_ptarcade",
        "SIGW-delta": "sigw_delta_ptarcade",
        "Cosmic superstrings": "super_ptarcade",
    }
    cosmology_rows = []
    for label, key in cosmology_map.items():
        agg = robust["aggregate"][key]
        cosmology_rows.append(
            {
                "label": label,
                "lnB_mean": agg["lnB_mean"],
                "sample_std": agg["lnB_sample_std"],
                "mean_nested_lnB_err": agg["mean_nested_lnB_err"],
                "published_lnB": agg["published_lnB"],
                "mean_residual_vs_published": agg["residual_mean"],
                "max_abs_residual_vs_published": agg["max_abs_residual"],
            }
        )

    h6_rows = h456["summary"]["H6_sigw_delta_prior"]
    sigw_delta_prior_values = [row["lnB"] for row in h6_rows]
    sigw_delta_prior_deltas = [
        row["delta_lnB_vs_original"]
        for row in h6_rows
        if row["delta_lnB_vs_original"] is not None
    ]

    qmc_rows = []
    for key in ["sigw_gauss", "sigw_delta", "super", "smbhb_env_fixed_gamma"]:
        item = qmc["aggregate"][key]
        qmc_rows.append(
            {
                "label": item["label"],
                "lnB_ti_mean": item.get("lnB_ti_mean"),
                "lnB_ti_sample_std": item.get("lnB_ti_sample_std"),
                "nested_robust_mean_lnB": item.get("nested_robust_mean_lnB"),
                "ti_residual_vs_nested_mean": item.get("ti_residual_vs_nested_mean"),
            }
        )

    h7_qmc = h7["qmc_ti_crosscheck"]["aggregate"]
    astro_rows = []
    for key, agg in sorted(
        h7["aggregate"].items(), key=lambda kv: kv[1]["lnB_mean"], reverse=True
    ):
        q = h7_qmc.get(key)
        astro_rows.append(
            {
                "model": key,
                "class": agg["class"],
                "lnB_mean": agg["lnB_mean"],
                "sample_std": agg["lnB_sample_std"],
                "mean_nested_lnB_err": agg["mean_nested_lnB_err"],
                "qmc_ti_lnB": None if q is None else q["lnB_ti_mean"],
                "qmc_ti_residual_vs_nested": None
                if q is None
                else float(q["lnB_ti_mean"] - agg["lnB_mean"]),
            }
        )

    class_best = {}
    for row in astro_rows:
        old = class_best.get(row["class"])
        if old is None or row["lnB_mean"] > old["lnB_mean"]:
            class_best[row["class"]] = row

    h9_summary = h9["summary_by_model"]
    bin_rows = []
    for key, item in sorted(
        h9_summary.items(), key=lambda kv: kv[1]["lnB_final_N14"], reverse=True
    ):
        bin_rows.append(
            {
                "model": key,
                "group": item["group"],
                "lnB_N14": item["lnB_final_N14"],
                "first_N_within_0p5nat_of_N14": item["first_N_within_0p5nat_of_N14"],
                "delta_lnB_N6_to_N14": item["delta_lnB_N6_to_N14"],
                "delta_lnB_N8_to_N14": item["delta_lnB_N8_to_N14"],
            }
        )

    cosmology_means = [row["lnB_mean"] for row in cosmology_rows]
    astro_class_best_means = [row["lnB_mean"] for row in class_best.values()]
    astro_all_means = [row["lnB_mean"] for row in astro_rows]
    top_cosmo = max(cosmology_rows, key=lambda row: row["lnB_mean"])
    top_astro = max(astro_rows, key=lambda row: row["lnB_mean"])

    overlap_low = max(min(cosmology_means), min(astro_class_best_means))
    overlap_high = min(max(cosmology_means), max(astro_class_best_means))
    overlap_exists = overlap_low <= overlap_high

    systematic_budget = {
        "official_density_reproduction": {
            "max_abs_mean_residual_vs_published": max_abs(
                [row["mean_residual_vs_published"] for row in cosmology_rows]
            ),
            "max_abs_single_config_residual_vs_published": max(
                row["max_abs_residual_vs_published"] for row in cosmology_rows
            ),
            "max_seed_livepoint_sample_std": max(row["sample_std"] for row in cosmology_rows),
        },
        "qmc_ti_crosscheck": {
            "diagnostic_pass": qmc["diagnostic_pass"],
            "max_abs_ti_residual_vs_dynesty_nested_mean": max_abs(
                [
                    row["ti_residual_vs_nested_mean"]
                    for row in qmc_rows
                    if row["ti_residual_vs_nested_mean"] is not None
                ]
            ),
            "max_abs_h7_ti_residual_vs_dynesty_nested_mean": max_abs(
                [
                    row["qmc_ti_residual_vs_nested"]
                    for row in astro_rows
                    if row["qmc_ti_residual_vs_nested"] is not None
                ]
            ),
        },
        "prior_boundary": {
            "sigw_delta_prior_lnb_min": min(sigw_delta_prior_values),
            "sigw_delta_prior_lnb_max": max(sigw_delta_prior_values),
            "sigw_delta_max_delta_lnb_vs_original": max(sigw_delta_prior_deltas),
        },
        "astrophysical_family": {
            "family_pass": h7["family_pass"],
            "family_pass_classes": h7["family_pass_classes"],
            "all_astro_lnB_min": min(astro_all_means),
            "all_astro_lnB_max": max(astro_all_means),
            "class_best_lnB_min": min(astro_class_best_means),
            "class_best_lnB_max": max(astro_class_best_means),
            "environmental_prior_variant_span": max(
                row["lnB_mean"]
                for row in astro_rows
                if row["class"] == "environmental_turnover"
            )
            - min(
                row["lnB_mean"]
                for row in astro_rows
                if row["class"] == "environmental_turnover"
            ),
        },
        "frequency_bins": {
            "models_within_0p5nat_by_N8": [
                row["model"]
                for row in bin_rows
                if row["first_N_within_0p5nat_of_N14"] <= 8
            ],
            "models_requiring_N10_for_0p5nat": [
                row["model"]
                for row in bin_rows
                if row["first_N_within_0p5nat_of_N14"] == 10
            ],
            "max_abs_delta_lnB_N8_to_N14": max_abs(
                [row["delta_lnB_N8_to_N14"] for row in bin_rows]
            ),
        },
        "source_identification": {
            "cosmology_robust_mean_range": [min(cosmology_means), max(cosmology_means)],
            "astrophysical_class_best_range": [
                min(astro_class_best_means),
                max(astro_class_best_means),
            ],
            "astro_cosmo_overlap_exists": overlap_exists,
            "astro_cosmo_overlap_range": [overlap_low, overlap_high]
            if overlap_exists
            else None,
            "top_cosmology": top_cosmo,
            "top_astrophysical_curvature": top_astro,
            "top_cosmo_minus_top_astro": float(
                top_cosmo["lnB_mean"] - top_astro["lnB_mean"]
            ),
            "sigw_delta_minus_top_astro": float(
                next(row for row in cosmology_rows if row["label"] == "SIGW-delta")[
                    "lnB_mean"
                ]
                - top_astro["lnB_mean"]
            ),
        },
    }

    return {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "task": "PRL H10 systematic evidence envelope",
        "source_files": {
            "official_density_robustness": str(ROBUST_JSON.relative_to(ROOT)),
            "h4_h5_h6": str(H456_JSON.relative_to(ROOT)),
            "qmc_ti": str(QMC_JSON.relative_to(ROOT)),
            "h7_astro_family": str(H7_JSON.relative_to(ROOT)),
            "h9_bin_driver": str(H9_JSON.relative_to(ROOT)),
        },
        "cosmology_rows": cosmology_rows,
        "sigw_delta_prior_rows": h6_rows,
        "qmc_ti_rows": qmc_rows,
        "astrophysical_family_rows": astro_rows,
        "astrophysical_class_best": class_best,
        "bin_driver_rows": bin_rows,
        "systematic_budget": systematic_budget,
    }


def write_md(data):
    budget = data["systematic_budget"]
    lines = [
        "# PRL Systematic Evidence Envelope",
        "",
        f"**Generated**: {data['generated']}",
        f"**JSON**: `{OUTFILE.relative_to(ROOT)}`",
        "",
        "## Cosmological Official-Density Reproduction",
        "",
        "| Model | mean lnB | sample std | mean residual | max abs residual |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in data["cosmology_rows"]:
        lines.append(
            f"| {row['label']} | `{row['lnB_mean']:+.3f}` | `{row['sample_std']:.3f}` | "
            f"`{row['mean_residual_vs_published']:+.3f}` | "
            f"`{row['max_abs_residual_vs_published']:.3f}` |"
        )
    lines.extend(
        [
            "",
            "## Curved SMBHB Control Family",
            "",
            "| Model | class | mean lnB | sample std | QMC/TI residual |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in data["astrophysical_family_rows"]:
        qres = row["qmc_ti_residual_vs_nested"]
        qtxt = "" if qres is None else f"`{qres:+.3f}`"
        lines.append(
            f"| `{row['model']}` | {row['class']} | `{row['lnB_mean']:+.3f}` | "
            f"`{row['sample_std']:.3f}` | {qtxt} |"
        )
    lines.extend(
        [
            "",
            "## Budget Summary",
            "",
            "| Component | Envelope / diagnostic |",
            "|---|---|",
            (
                "| Official-density reproduction | "
                f"mean residual max `{budget['official_density_reproduction']['max_abs_mean_residual_vs_published']:.3f}` nat; "
                f"single-config max `{budget['official_density_reproduction']['max_abs_single_config_residual_vs_published']:.3f}` nat |"
            ),
            (
                "| QMC/TI cross-check | "
                f"main max residual `{budget['qmc_ti_crosscheck']['max_abs_ti_residual_vs_dynesty_nested_mean']:.3f}` nat; "
                f"curved-family max residual `{budget['qmc_ti_crosscheck']['max_abs_h7_ti_residual_vs_dynesty_nested_mean']:.3f}` nat |"
            ),
            (
                "| SIGW-delta prior boundary | "
                f"`Delta lnB <= {budget['prior_boundary']['sigw_delta_max_delta_lnb_vs_original']:.3f}` nat over scanned upper bounds |"
            ),
            (
                "| Curved SMBHB family | "
                f"tested-representative range `{budget['astrophysical_family']['class_best_lnB_min']:+.3f}` to "
                f"`{budget['astrophysical_family']['class_best_lnB_max']:+.3f}`; "
                f"family pass `{budget['astrophysical_family']['family_pass']}` |"
            ),
            (
                "| Frequency-bin driver | "
                f"max `|Delta lnB(N=8 to 14)| = {budget['frequency_bins']['max_abs_delta_lnB_N8_to_N14']:.3f}` nat; "
                f"most models within 0.5 nat by `N=8` |"
            ),
            "",
            "## Source-Identification Consequence",
            "",
        ]
    )
    sid = budget["source_identification"]
    lines.append(
        "The robust cosmological means span "
        f"`{sid['cosmology_robust_mean_range'][0]:+.3f}` to "
        f"`{sid['cosmology_robust_mean_range'][1]:+.3f}`.  The best row in each "
        "curved-SMBHB class spans "
        f"`{sid['astrophysical_class_best_range'][0]:+.3f}` to "
        f"`{sid['astrophysical_class_best_range'][1]:+.3f}`.  "
        f"The ranges overlap: `{sid['astro_cosmo_overlap_exists']}`."
    )
    lines.append("")
    lines.append(
        "The top cosmological mean exceeds the top curved-SMBHB mean by "
        f"`{sid['top_cosmo_minus_top_astro']:.3f}` nat, while SIGW-delta exceeds the "
        f"top curved-SMBHB row by only `{sid['sigw_delta_minus_top_astro']:.3f}` nat. "
        "Therefore the PRL claim should be framed as non-discrimination on the "
        "current calibrated evidence scale, not as a new-physics source detection."
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main():
    data = build_envelope()
    OUTFILE.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    write_md(data)
    print(f"saved: {OUTFILE}")
    print(f"saved: {OUT_MD}")


if __name__ == "__main__":
    main()
