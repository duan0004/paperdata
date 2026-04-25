#!/usr/bin/env python3
"""Final PRL stability diagnostics assembled from completed outputs.

This script does not run new evidence calculations.  It consolidates existing
official-density, bridge, PPC, and SMBHB-control outputs into small tables that
support the final PRL framing: current public PTA data identify low-frequency
curvature more robustly than a unique source class.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results" / "T2_NG15yr" / "bayes_factors"
OUT_JSON = OUTDIR / "prl_final_stability_diagnostics.json"
OUT_MD = OUTDIR / "prl_final_stability_diagnostics.md"

H10 = OUTDIR / "prl_H10_systematic_envelope.json"
IDENT = ROOT / "results" / "identifiability"
PPC = ROOT / "results" / "ppc"
PHYS = ROOT / "results" / "physical_smbhb" / "ng15_official_family.csv"
ROBUST = ROOT / "results" / "prl_reference_bridge" / "robustness_budget.csv"


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


def f3(x: float) -> str:
    return f"{float(x):.3f}"


def bridge_seed_livepoint_range() -> float:
    rows = read_csv(ROBUST)
    by_model: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        if row.get("status") != "pass":
            continue
        by_model[row["model"]].append(float(row["ln_Z"]))
    ranges = [max(vals) - min(vals) for vals in by_model.values() if vals]
    return max(ranges) if ranges else float("nan")


def family_ppc_summary() -> list[dict]:
    low_rows = read_csv(PPC / "family_ppc_lowfreq.csv")
    bin_rows = read_csv(PPC / "posterior_summary_residuals.csv")
    out = []
    for family in sorted({r["family"] for r in low_rows}):
        low_vals = [float(r["mean_abs_z_low4"]) for r in low_rows if r["family"] == family]
        cov_vals = [float(r["coverage_low4"]) for r in low_rows if r["family"] == family]
        max_low = max(float(r["max_abs_z_low4"]) for r in low_rows if r["family"] == family)
        high_vals = [
            abs(float(r["z_median"]))
            for r in bin_rows
            if r["family"] == family and int(r["bin"]) >= 9
        ]
        out.append(
            {
                "family": family,
                "mean_abs_z_low4": sum(low_vals) / len(low_vals),
                "max_abs_z_low4": max_low,
                "coverage_low4": sum(cov_vals) / len(cov_vals),
                "mean_abs_z_high_bins": sum(high_vals) / len(high_vals),
                "max_abs_z_high_bins": max(high_vals),
                "interpretation": "posterior-summary diagnostic only; comparable low-frequency residuals across top families",
            }
        )
    return out


def environmental_prior_rows(top_cosmo: float, sigw_delta: float, superstrings: float) -> list[dict]:
    rows = []
    for row in read_csv(PHYS):
        if row["class"] != "environmental_turnover":
            continue
        lnb = float(row["mean_lnB_vs_ptarcade_bhb"])
        rows.append(
            {
                "model": row["model_key"],
                "physical_status": row["physical_status"],
                "lnB": lnb,
                "sample_std": float(row["sample_std"]),
                "gap_to_sigw_gaussian": top_cosmo - lnb,
                "gap_to_sigw_delta": sigw_delta - lnb,
                "gap_to_superstrings": superstrings - lnb,
                "interpretation": row["claim_scope"],
            }
        )
    return sorted(rows, key=lambda r: r["lnB"], reverse=True)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    h10 = json.loads(H10.read_text())
    budget = h10["systematic_budget"]
    source_id = budget["source_identification"]
    top_cosmo = float(source_id["top_cosmology"]["lnB_mean"])
    top_astro = float(source_id["top_astrophysical_curvature"]["lnB_mean"])
    sigw_delta_gap = float(source_id["sigw_delta_minus_top_astro"])
    top_gap = float(source_id["top_cosmo_minus_top_astro"])
    sigw_delta = top_astro + sigw_delta_gap
    superstrings = float(source_id["cosmology_robust_mean_range"][0])

    family_gaps = read_csv(IDENT / "family_gap_summary.csv")
    basis = read_csv(IDENT / "spectral_basis_evidence.csv")
    gamma = read_csv(IDENT / "gamma_projection_summary.csv")
    ppc = family_ppc_summary()
    env_rows = environmental_prior_rows(top_cosmo, sigw_delta, superstrings)
    bridge_range = bridge_seed_livepoint_range()

    astro_family_delta = next(
        float(r["delta_to_best_family"]) for r in family_gaps if r["family"] == "Astro-curved"
    )
    gamma_min = min(float(r["gamma_eff_min_median"]) for r in gamma)
    gamma_max = max(float(r["gamma_eff_max_median"]) for r in gamma)

    uncertainty_rows = [
        {
            "effect": "official-density reproduction residual",
            "size_nat": budget["official_density_reproduction"]["max_abs_single_config_residual_vs_published"],
            "scope": "maximum absolute single-configuration residual among the three calibrated NANOGrav targets",
        },
        {
            "effect": "SIGW-delta peak-frequency prior boundary",
            "size_nat": budget["prior_boundary"]["sigw_delta_max_delta_lnb_vs_original"],
            "scope": "scan over tested upper peak-frequency bounds",
        },
        {
            "effect": "dynesty/QMC-TI residual, main rows",
            "size_nat": budget["qmc_ti_crosscheck"]["max_abs_ti_residual_vs_dynesty_nested_mean"],
            "scope": "independent Sobol-QMC/TI cross-check",
        },
        {
            "effect": "dynesty/QMC-TI residual, curved rows",
            "size_nat": budget["qmc_ti_crosscheck"]["max_abs_h7_ti_residual_vs_dynesty_nested_mean"],
            "scope": "family cross-check for tested curved-SMBHB representatives",
        },
        {
            "effect": "hybrid3 leading-row seed/live range",
            "size_nat": bridge_range,
            "scope": "maximum range over the bridge robustness sweep",
        },
        {
            "effect": "SIGW-Gaussian minus best curved-SMBHB model",
            "size_nat": top_gap,
            "scope": "matched NG15 official-density scale",
        },
        {
            "effect": "SIGW-like minus Astro-curved family",
            "size_nat": abs(astro_family_delta),
            "scope": "hybrid3 equal-weight tested-representative family mixture",
        },
    ]

    data = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "source_files": {
            "h10": str(H10.relative_to(ROOT)),
            "identifiability": str((IDENT / "mode_compression.md").relative_to(ROOT)),
            "ppc": str((PPC / "family_ppc_lowfreq.csv").relative_to(ROOT)),
            "physical_smbhb": str(PHYS.relative_to(ROOT)),
            "bridge_robustness": str(ROBUST.relative_to(ROOT)),
        },
        "mode_compression": {
            "spectral_basis_rows": basis,
            "astro_curved_family_delta_nat": astro_family_delta,
            "gamma_eff_median_range": [gamma_min, gamma_max],
            "interpretation": "amplitude-like and low-frequency-curvature diagnostics are more stable than source-family labels",
        },
        "family_ppc": ppc,
        "uncertainty_budget": uncertainty_rows,
        "environmental_prior_sensitivity": env_rows,
    }
    OUT_JSON.write_text(json.dumps(data, indent=2) + "\n")

    write_csv(
        OUTDIR / "prl_final_family_ppc_summary.csv",
        ppc,
        [
            "family",
            "mean_abs_z_low4",
            "max_abs_z_low4",
            "coverage_low4",
            "mean_abs_z_high_bins",
            "max_abs_z_high_bins",
            "interpretation",
        ],
    )
    write_csv(
        OUTDIR / "prl_final_uncertainty_budget.csv",
        uncertainty_rows,
        ["effect", "size_nat", "scope"],
    )
    write_csv(
        OUTDIR / "prl_final_environmental_prior_sensitivity.csv",
        env_rows,
        [
            "model",
            "physical_status",
            "lnB",
            "sample_std",
            "gap_to_sigw_gaussian",
            "gap_to_sigw_delta",
            "gap_to_superstrings",
            "interpretation",
        ],
    )

    lines = [
        "# PRL Final Stability Diagnostics",
        "",
        f"Generated: {data['generated']}",
        "",
        "These tables consolidate existing completed runs. They do not introduce new timing-level or nested-sampling evidence.",
        "",
        "## Curvature-Mode Compression",
        "",
        "| diagnostic | value | interpretation |",
        "|---|---:|---|",
        f"| Astro-curved family gap on hybrid3 | `{astro_family_delta:+.3f}` nat | source-family separation remains order-unity |",
        f"| curved-SMBHB gamma_eff median range | `{gamma_min:.3f}`--`{gamma_max:.3f}` | window-dependent low-frequency curvature is directly visible |",
        f"| top model gap, NG15 official scale | `{top_gap:.3f}` nat | comparable to calibration/prior systematics |",
        "",
        "## Posterior-Summary PPC",
        "",
        "| family | mean low-4 | max low-4 | mean high-bin | interpretation |",
        "|---|---:|---:|---:|---|",
    ]
    for row in ppc:
        lines.append(
            f"| {row['family']} | `{f3(row['mean_abs_z_low4'])}` | `{f3(row['max_abs_z_low4'])}` | "
            f"`{f3(row['mean_abs_z_high_bins'])}` | {row['interpretation']} |"
        )
    lines.extend(
        [
            "",
            "## Evidence/Systematics Budget",
            "",
            "| effect | size | scope |",
            "|---|---:|---|",
        ]
    )
    for row in uncertainty_rows:
        lines.append(f"| {row['effect']} | `{f3(row['size_nat'])}` nat | {row['scope']} |")
    lines.extend(
        [
            "",
            "## Environmental Prior Sensitivity",
            "",
            "| model | lnB | gap to SIGW-Gaussian | gap to SIGW-delta | gap to superstrings | status |",
            "|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in env_rows:
        lines.append(
            f"| `{row['model']}` | `{f3(row['lnB'])}` | `{f3(row['gap_to_sigw_gaussian'])}` | "
            f"`{f3(row['gap_to_sigw_delta'])}` | `{f3(row['gap_to_superstrings'])}` | {row['physical_status']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The final stability diagnostics support the PRL framing without upgrading the claim to full timing-level source identification. The public posterior-summary PPCs show comparable low-frequency residuals for SIGW-like, Astro-curved, and Cosmo-other leading families; the matched-scale evidence gaps remain comparable to identified calibration, prior-boundary, and bridge robustness systematics. Environmental-turnover prior variants remain competitive with at least part of the leading cosmological tier, but they do not constitute a complete population-synthesis marginalization.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")
    print({"json": str(OUT_JSON), "md": str(OUT_MD), "ppc_families": len(ppc)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
