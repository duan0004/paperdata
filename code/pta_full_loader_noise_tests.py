#!/usr/bin/env python3
"""Full public-PTA loader and noise-product tests.

This is a read-only production gate.  It uses ENTERPRISE/enterprise_extensions
loaders, never edits `.par` timing models, and records array-specific noise
products that are present locally.  The output is a loader/noise integration
manifest, not a timing-level Bayes-factor calculation.
"""

from __future__ import annotations

import csv
import json
import sys
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
OUT_DIR = ROOT / "results" / "5pta_timing"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

from pta_loader_smoke import (  # noqa: E402
    epta_pairs,
    inpta_pairs,
    load_epta,
    load_generic_par_tim,
    load_ng15,
    load_ppta,
    mpta_pairs,
    ppta_pairs,
    rel,
)


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def load_json_len(path: Path) -> int | str:
    if not path.exists():
        return "missing"
    try:
        with path.open() as f:
            obj = json.load(f)
        return len(obj) if hasattr(obj, "__len__") else "scalar"
    except Exception as exc:
        return f"{type(exc).__name__}: {exc}"


def noise_inventory() -> list[dict]:
    rows = []
    ng_noise = ROOT / "data" / "NG15yr" / "tutorials" / "data" / "15yr_wn_dict.json"
    rows.append(
        {
            "array": "NG15",
            "noise_root": rel(ng_noise),
            "noise_json_files": int(ng_noise.exists()),
            "noise_param_entries": load_json_len(ng_noise),
            "auxiliary_entries": "",
            "note": "Official NANOGrav 15yr white-noise dictionary used by local NG15 loaders.",
        }
    )

    epta_root = ROOT / "data" / "EPTA_DR2" / "epta-dr2" / "EPTA-DR2" / "noisefiles" / "DR2new+"
    rows.append(
        {
            "array": "EPTA_DR2new+",
            "noise_root": rel(epta_root),
            "noise_json_files": len(list(epta_root.glob("*_noise.json"))) if epta_root.exists() else 0,
            "noise_param_entries": "",
            "auxiliary_entries": len(list(epta_root.glob("*_dict.json"))) if epta_root.exists() else 0,
            "note": "EPTA DR2new+ per-pulsar noise JSON plus red/dm/chrom dictionaries.",
        }
    )

    ppta_root = ROOT / "data" / "PPTA_DR3" / "ppta_dr3" / "toas_and_parameters" / "noisefiles"
    rows.append(
        {
            "array": "PPTA_DR3",
            "noise_root": rel(ppta_root),
            "noise_json_files": len(list(ppta_root.glob("*_noise.json"))) if ppta_root.exists() else 0,
            "noise_param_entries": "",
            "auxiliary_entries": "",
            "note": "PPTA DR3 single-pulsar noise JSON files.",
        }
    )

    inpta_root = ROOT / "data" / "InPTA_DR1" / "InPTA.DR1"
    rows.append(
        {
            "array": "InPTA_DR1",
            "noise_root": rel(inpta_root),
            "noise_json_files": 0,
            "noise_param_entries": "",
            "auxiliary_entries": len(list(inpta_root.glob("*/*DM_timeseries*"))) if inpta_root.exists() else 0,
            "note": "No explicit noise JSON was located in the staged public repository; DM time-series products are present.",
        }
    )

    refreshed_mpta_root = ROOT / "data" / "MPTA_4p5yr" / "partim_datacentral_20260425" / "partim"
    legacy_mpta_root = ROOT / "data" / "MPTA_4p5yr" / "partim" / "noise_subtracted_residuals" / "DOI_partim"
    mpta_root = refreshed_mpta_root if refreshed_mpta_root.exists() else legacy_mpta_root
    rows.append(
        {
            "array": "MPTA_4p5yr",
            "noise_root": rel(mpta_root),
            "noise_json_files": 0,
            "noise_param_entries": "",
            "auxiliary_entries": "",
            "note": "Staged public product is the noise-subtracted residual par/tim release; no separate noise JSON was located.",
        }
    )
    return rows


def status_for(res: dict) -> str:
    if res.get("status") == "BLOCKED":
        return "BLOCKED"
    available = int(res.get("available_pulsars", 0))
    loaded = len(res.get("loaded_pulsars", []))
    failures = len(res.get("failures", []))
    if available > 0 and loaded == available and failures == 0:
        return "PASS"
    if loaded > 0 and failures > 0:
        return "PARTIAL"
    return "FAIL"


def run_loaders() -> list[dict]:
    configs = [
        ("NG15", lambda: load_ng15(0), len(load_ng15(0).get("loaded_pulsars", []))),
        ("EPTA_DR2new+", lambda: load_epta(0, "DR2new+", "FROM_PAR", None), len(epta_pairs("DR2new+"))),
        ("PPTA_DR3", lambda: load_ppta(0, "FROM_PAR", None), len(ppta_pairs())),
        (
            "InPTA_DR1_NB",
            lambda: load_generic_par_tim("InPTA_DR1_NB", inpta_pairs("NB"), 0, "FROM_PAR", None),
            len(inpta_pairs("NB")),
        ),
        (
            "InPTA_DR1_WB",
            lambda: load_generic_par_tim("InPTA_DR1_WB", inpta_pairs("WB"), 0, "FROM_PAR", None),
            len(inpta_pairs("WB")),
        ),
        ("MPTA_4p5yr", lambda: load_generic_par_tim("MPTA_4p5yr", mpta_pairs(), 0, "FROM_PAR", None), len(mpta_pairs())),
    ]
    rows = []
    for label, fn, expected_pairs in configs:
        try:
            res = fn()
        except Exception as exc:
            res = {
                "array": label,
                "status": "FAIL",
                "available_pulsars": expected_pairs,
                "loaded_pulsars": [],
                "failures": [{"pulsar": "", "traceback": f"{type(exc).__name__}: {exc}"}],
            }
        res["test_label"] = label
        res["expected_pairs"] = expected_pairs
        res["status"] = status_for(res)
        path = OUT_DIR / f"full_loader_noise_{label}.json"
        path.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
        rows.append(
            {
                "test_label": label,
                "array": res.get("array", label),
                "status": res["status"],
                "expected_pairs": expected_pairs,
                "available_pulsars": res.get("available_pulsars", ""),
                "loaded": len(res.get("loaded_pulsars", [])),
                "failures": len(res.get("failures", [])),
                "json": rel(path),
            }
        )
        print(f"{label}: {res['status']} loaded={rows[-1]['loaded']} failures={rows[-1]['failures']}", flush=True)
    return rows


def write_markdown(loader_rows: list[dict], noise_rows: list[dict]) -> None:
    lines = [
        "# Full Public-PTA Loader And Noise Tests",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Scope",
        "",
        "- Read-only ENTERPRISE/enterprise_extensions loader gate.",
        "- No timing model `.par` files are edited.",
        "- This is not a timing-level evidence calculation.",
        "",
        "## Loader Results",
        "",
        "| test | status | expected pairs | loaded | failures | JSON |",
        "|---|---|---:|---:|---:|---|",
    ]
    for row in loader_rows:
        lines.append(
            f"| `{row['test_label']}` | `{row['status']}` | {row['expected_pairs']} | "
            f"{row['loaded']} | {row['failures']} | `{row['json']}` |"
        )
    lines.extend(
        [
            "",
            "## Noise Product Inventory",
            "",
            "| array | noise root | noise JSON files | auxiliary entries | note |",
            "|---|---|---:|---:|---|",
        ]
    )
    for row in noise_rows:
        lines.append(
            f"| `{row['array']}` | `{row['noise_root']}` | {row['noise_json_files']} | "
            f"{row['auxiliary_entries']} | {row['note']} |"
        )
    lines.extend(
        [
            "",
            "## Gate Interpretation",
            "",
            "- PASS means every discovered par/tim or feather product in that test loaded through the standard loader.",
            "- PARTIAL means the public product is staged but at least one loader edge case must be resolved or excluded before production evidence.",
            "- Arrays without explicit noise JSON require a documented likelihood policy before a full five-PTA evidence claim.",
            "",
        ]
    )
    (OUT_DIR / "full_loader_noise_tests.md").write_text("\n".join(lines))


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    loader_rows = run_loaders()
    noise_rows = noise_inventory()
    write_csv(
        OUT_DIR / "full_loader_noise_summary.csv",
        loader_rows,
        ["test_label", "array", "status", "expected_pairs", "available_pulsars", "loaded", "failures", "json"],
    )
    write_csv(
        OUT_DIR / "noise_product_inventory.csv",
        noise_rows,
        ["array", "noise_root", "noise_json_files", "noise_param_entries", "auxiliary_entries", "note"],
    )
    write_markdown(loader_rows, noise_rows)
    print({"loader_tests": len(loader_rows), "noise_rows": len(noise_rows), "report": rel(OUT_DIR / "full_loader_noise_tests.md")})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
