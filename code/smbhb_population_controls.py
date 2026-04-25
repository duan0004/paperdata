#!/usr/bin/env python3
"""P5 physical-status audit for SMBHB control models.

This consolidates existing official-density and bridge results into a
population-control status table.  It does not pretend that phenomenological
curvature controls are population-synthesis models.
"""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "physical_smbhb"
THEORY = ROOT / "theory" / "SMBHB_population_controls.md"
H7 = ROOT / "results" / "T2_NG15yr" / "bayes_factors" / "prl_H7_astrophysical_family.json"
TIMING = ROOT / "results" / "5pta_timing" / "template_lnZ.csv"


PHYSICAL_STATUS = {
    "env_fixed_gamma": {
        "label": "environmental turnover, fixed gamma",
        "status": "physical control",
        "claim_scope": "SMBHB environmental-hardening spectral control",
        "reference": "Sampson et al. 2015, arXiv:1505.07179; Sesana 2013, arXiv:1211.5375",
    },
    "env_free_gamma": {
        "label": "environmental turnover, free gamma",
        "status": "physical control with extra spectral freedom",
        "claim_scope": "environmental-hardening control plus free effective gamma",
        "reference": "Sampson et al. 2015, arXiv:1505.07179",
    },
    "env_lowbend_fixed_gamma": {
        "label": "environmental turnover, low-bend prior",
        "status": "physical-prior sensitivity",
        "claim_scope": "environmental-hardening prior-range sensitivity",
        "reference": "Sampson et al. 2015, arXiv:1505.07179",
    },
    "env_broadbend_fixed_gamma": {
        "label": "environmental turnover, broad-bend prior",
        "status": "physical-prior sensitivity",
        "claim_scope": "environmental-hardening prior-range sensitivity",
        "reference": "Sampson et al. 2015, arXiv:1505.07179",
    },
    "broken_pl_fixed_gamma": {
        "label": "broken-power-law curvature",
        "status": "phenomenological curvature control",
        "claim_scope": "low-frequency suppression surrogate, not population synthesis",
        "reference": "Sampson et al. 2015, arXiv:1505.07179",
    },
    "ecc_supp_fixed_gamma": {
        "label": "eccentricity-inspired suppression",
        "status": "phenomenological eccentricity-inspired control",
        "claim_scope": "effective low-frequency suppression, not population synthesis",
        "reference": "Sampson et al. 2015, arXiv:1505.07179",
    },
}


def read_csv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    h7 = json.loads(H7.read_text())
    aggregate = h7["aggregate"]
    qmc = h7.get("qmc_ti_crosscheck", {})
    qmc_aggregate = qmc.get("aggregate", {})
    rows = []
    for key, info in PHYSICAL_STATUS.items():
        if key not in aggregate:
            continue
        agg = aggregate[key]
        q = qmc_aggregate.get(key, {})
        rows.append(
            {
                "model_key": key,
                "label": info["label"],
                "class": agg.get("class", ""),
                "physical_status": info["status"],
                "claim_scope": info["claim_scope"],
                "mean_lnB_vs_ptarcade_bhb": f"{float(agg['lnB_mean']):.12g}",
                "sample_std": f"{float(agg['lnB_sample_std']):.12g}",
                "min_lnB": f"{float(agg['lnB_min']):.12g}",
                "max_lnB": f"{float(agg['lnB_max']):.12g}",
                "mean_nested_err": f"{float(agg['mean_nested_lnB_err']):.12g}",
                "qmc_ti_lnB": "" if not q else f"{float(q['lnB_ti_mean']):.12g}",
                "reference": info["reference"],
            }
        )
    write_csv(
        OUT / "ng15_official_family.csv",
        rows,
        [
            "model_key",
            "label",
            "class",
            "physical_status",
            "claim_scope",
            "mean_lnB_vs_ptarcade_bhb",
            "sample_std",
            "min_lnB",
            "max_lnB",
            "mean_nested_err",
            "qmc_ti_lnB",
            "reference",
        ],
    )

    timing_rows = []
    if TIMING.exists():
        for r in read_csv(TIMING):
            if r["family"] not in {"Astro-curved", "Astro-simple"}:
                continue
            timing_rows.append(
                {
                    "model": r["model"],
                    "family": r["family"],
                    "evidence_scale": r["evidence_scale"],
                    "full_timing_bayes_factor": r["full_timing_bayes_factor"],
                    "ln_Z": r["ln_Z"],
                    "ln_Z_err": r["ln_Z_err"],
                    "delta_to_best": r["delta_to_best"],
                    "interpretation": "bridge stress-test row; not full timing-level Bayes factor",
                }
            )
    write_csv(
        OUT / "timing_family.csv",
        timing_rows,
        [
            "model",
            "family",
            "evidence_scale",
            "full_timing_bayes_factor",
            "ln_Z",
            "ln_Z_err",
            "delta_to_best",
            "interpretation",
        ],
    )

    lines = [
        "# SMBHB Population-Control Status",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Scope",
        "",
        "P5 separates physically interpretable SMBHB controls from phenomenological curvature controls. The current PRL evidence supports a calibrated curved-SMBHB control statement, not a complete SMBHB population-synthesis marginalization.",
        "",
        "## NG15 Official-Density Scale",
        "",
        "| Model | physical status | lnB vs PTArcade BHB | scope |",
        "|---|---|---:|---|",
    ]
    for r in rows:
        lines.append(
            f"| `{r['label']}` | `{r['physical_status']}` | {float(r['mean_lnB_vs_ptarcade_bhb']):+.3f} +/- {float(r['sample_std']):.3f} | {r['claim_scope']} |"
        )
    lines.extend(
        [
            "",
            "## Bridge/Timing-Side Status",
            "",
            "| Model | family | ln Z | delta to best | status |",
            "|---|---|---:|---:|---|",
        ]
    )
    for r in timing_rows:
        lines.append(
            f"| `{r['model']}` | `{r['family']}` | {float(r['ln_Z']):+.3f} +/- {float(r['ln_Z_err']):.3f} | {float(r['delta_to_best']):+.3f} | {r['interpretation']} |"
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Main-text safe: environmental turnover as a physical SMBHB control; curved-SMBHB family as a tested-representative control family.",
            "- Supplement safe: broken-power-law and eccentricity-inspired controls as low-frequency curvature surrogates.",
            "- Not yet safe: claiming full astrophysical population-synthesis marginalization. That requires external population priors or simulation-derived template weights and a rerun through the official-density/timing-level evidence gates.",
        ]
    )
    THEORY.write_text("\n".join(lines) + "\n")
    print({"ng15_rows": len(rows), "timing_rows": len(timing_rows), "theory": str(THEORY)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
