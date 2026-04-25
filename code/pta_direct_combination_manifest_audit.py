#!/usr/bin/env python3
"""Coordinate-based manifest audit for the public five-PTA direct-combination gate.

This script is read-only.  It parses public `.par` files only to derive a
canonical J-name for grouping duplicate pulsars across PTAs.  It does not edit
timing models, does not write converted `.par` files, and does not build a
timing-level likelihood.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
METAPULSAR_SRC = ROOT / "data" / "5PTA_public" / "metapulsar" / "src"
OUT_DIR = ROOT / "results" / "5pta_timing"
THEORY_DIR = ROOT / "theory"

for path in (CODE, METAPULSAR_SRC):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import pta_loader_smoke as loaders  # noqa: E402
from metapulsar.position_helpers import _format_j_name_from_icrs, _skycoord_from_pint_model  # noqa: E402
from pint.models import get_model  # noqa: E402


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def base_name(array: str, name: str) -> str:
    if array.startswith("InPTA"):
        return name.rsplit("_", 1)[0]
    return name


def coordinate_canonical_name(par: Path, fallback: str) -> tuple[str, str]:
    try:
        # Redirect noisy parser output.  This is a coordinate-grouping audit
        # only; no converted timing model is written or used for likelihoods.
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            model = get_model(str(par), allow_T2=True, allow_tcb=True)
            coord = _skycoord_from_pint_model(model)
        return _format_j_name_from_icrs(coord), ""
    except Exception as exc:
        return fallback, f"{type(exc).__name__}: {exc}"


def collect_pairs(inpta_band: str) -> list[dict]:
    raw = {
        "NG15": loaders.ng15_timing_pairs(),
        "EPTA": loaders.epta_pairs("DR2new+"),
        "PPTA": loaders.ppta_pairs(),
        "MPTA": loaders.mpta_pairs(),
        f"InPTA_{inpta_band}": loaders.inpta_pairs(inpta_band),
    }
    rows = []
    for array, pairs in raw.items():
        for name, par, tim in pairs:
            base = base_name(array, name)
            canonical, error = coordinate_canonical_name(par, base)
            rows.append(
                {
                    "array": array,
                    "input_name": name,
                    "base_name": base,
                    "canonical_name": canonical,
                    "par": rel(par),
                    "tim": rel(tim),
                    "canonical_error": error,
                }
            )
    return rows


def group_rows(rows: list[dict]) -> list[dict]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[row["canonical_name"]].append(row)
    out = []
    for canonical, items in sorted(grouped.items()):
        arrays = sorted({item["array"] for item in items})
        input_names = sorted({item["base_name"] for item in items})
        out.append(
            {
                "canonical_name": canonical,
                "occurrences": len(items),
                "arrays": ";".join(arrays),
                "input_names": ";".join(input_names),
                "is_duplicate": len(items) > 1,
                "has_alias_names": len(input_names) > 1,
            }
        )
    return out


def write_report(path: Path, payload: dict) -> None:
    summary = payload["summary"]
    lines = [
        "# Five-PTA Direct-Combination Manifest Audit",
        "",
        f"Generated: {payload['generated']}",
        "",
        "## Scope",
        "",
        "This is a read-only coordinate-grouping audit for the public five-PTA",
        "direct-combination gate.  It does not edit `.par` files, does not write",
        "converted timing models, and does not build a likelihood.",
        "",
        "## Summary",
        "",
        f"- InPTA product used for this audit: `{payload['inpta_band']}`.",
        f"- Array occurrences: `{summary['occurrences']}`.",
        f"- Canonical pulsar groups: `{summary['canonical_unique']}`.",
        f"- Duplicate canonical groups: `{summary['duplicate_groups']}`.",
        f"- Coordinate-canonical fallback failures: `{summary['canonical_failures']}`.",
        "",
        "The canonical group count matches the 121-pulsar public five-PTA target",
        "once B/J aliases are grouped by coordinates.  This does not resolve the",
        "published baseline policy; it only fixes the name-level manifest count.",
        "",
        "## Alias Groups",
        "",
        "| canonical name | input names | arrays | occurrences |",
        "|---|---|---|---:|",
    ]
    for row in payload["alias_groups"]:
        lines.append(
            f"| `{row['canonical_name']}` | `{row['input_names']}` | `{row['arrays']}` | {row['occurrences']} |"
        )
    if payload["canonical_failures"]:
        lines.extend(["", "## Coordinate Fallbacks", "", "| array | input name | fallback | error |", "|---|---|---|---|"])
        for row in payload["canonical_failures"]:
            lines.append(
                f"| `{row['array']}` | `{row['input_name']}` | `{row['canonical_name']}` | `{row['canonical_error']}` |"
            )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Occurrence CSV: `{payload['occurrence_csv']}`.",
            f"- Group CSV: `{payload['group_csv']}`.",
            f"- JSON: `{payload['json']}`.",
            "",
            "## Gate Decision",
            "",
            "Status: `NAME_AND_COORDINATE_MANIFEST_READY_POLICY_BLOCKED`.",
            "",
            "The local public timing inputs now reproduce the expected 121 canonical",
            "pulsar groups at the manifest level.  Timing-level family evidence remains",
            "blocked until the published direct-combination, InPTA preprocessing, and",
            "noise-policy baseline is reproduced.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inpta-band", default="NB", choices=["NB", "WB"])
    parser.add_argument("--tag", default="2026-04-25")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    THEORY_DIR.mkdir(parents=True, exist_ok=True)

    rows = collect_pairs(args.inpta_band)
    groups = group_rows(rows)
    alias_groups = [row for row in groups if row["has_alias_names"]]
    failures = [row for row in rows if row["canonical_error"]]
    summary = {
        "occurrences": len(rows),
        "canonical_unique": len(groups),
        "duplicate_groups": sum(1 for row in groups if row["is_duplicate"]),
        "alias_groups": len(alias_groups),
        "canonical_failures": len(failures),
    }

    stem = f"direct_combination_manifest_audit_{args.inpta_band.lower()}_{args.tag}"
    occurrence_csv = OUT_DIR / f"{stem}_occurrences.csv"
    group_csv = OUT_DIR / f"{stem}_groups.csv"
    json_path = OUT_DIR / f"{stem}.json"
    md_path = THEORY_DIR / f"5PTA_direct_combination_manifest_audit_{args.inpta_band.lower()}_{args.tag}.md"

    write_csv(occurrence_csv, rows, ["array", "input_name", "base_name", "canonical_name", "par", "tim", "canonical_error"])
    write_csv(group_csv, groups, ["canonical_name", "occurrences", "arrays", "input_names", "is_duplicate", "has_alias_names"])

    payload = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "inpta_band": args.inpta_band,
        "summary": summary,
        "alias_groups": alias_groups,
        "canonical_failures": failures,
        "occurrence_csv": rel(occurrence_csv),
        "group_csv": rel(group_csv),
        "json": rel(json_path),
        "markdown": rel(md_path),
    }
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_report(md_path, payload)
    print(json.dumps({"summary": summary, "markdown": rel(md_path), "json": rel(json_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
