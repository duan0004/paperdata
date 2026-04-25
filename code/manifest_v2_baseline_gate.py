#!/usr/bin/env python3
"""Manifest-v2 and exact-baseline gate for public PTA timing work.

This script inventories local public timing inputs and writes the next-stage
gate for exact timing-level baseline reproduction.  It does not run a
likelihood and does not produce source-family evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tarfile
import zipfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = ROOT / "code"
OUT_DIR = ROOT / "results" / "5pta_timing"
THEORY_DIR = ROOT / "theory"

if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import pta_loader_smoke as loaders  # noqa: E402


REFERENCE_BASELINE = {
    "source": "arXiv:2512.08666 source lines 338-350, 359-428",
    "dataset_counts": {
        "EPTA": 25,
        "InPTA_DR1": 14,
        "MPTA_4p5yr": 83,
        "NG15": 68,
        "PPTA": 32,
    },
    "direct_combination_priority": "NG15 > EPTA > PPTA > MPTA; InPTA omitted from priority order because its pulsars are observed by other PTAs in the reference analysis.",
    "noise_policy": "NG15/PPTA/MPTA and seven EPTA backends use EFAC/EQUAD/ECORR; InPTA and the remaining EPTA backends use EFAC/EQUAD.",
}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def strip_band(name: str) -> str:
    for suffix in ("_NB", "_WB"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def ng15_names() -> list[dict]:
    import contextlib
    import io

    from enterprise_extensions.load_feathers import load_feathers_from_folder

    feather_dir = ROOT / "data" / "NG15yr" / "tutorials" / "data" / "feathers"
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        psrs = load_feathers_from_folder(str(feather_dir))
    return [
        {"array": "NG15", "name": p.name, "base_name": p.name, "source": "feather", "parfile": "", "timfile": ""}
        for p in psrs
    ]


def par_tim_rows(array: str, pairs: list[tuple[str, Path, Path]]) -> list[dict]:
    rows = []
    for name, par, tim in pairs:
        rows.append(
            {
                "array": array,
                "name": name,
                "base_name": strip_band(name),
                "source": "par_tim",
                "parfile": rel(par),
                "timfile": rel(tim),
            }
        )
    return rows


def local_inventory(epta_variant: str) -> dict[str, list[dict]]:
    return {
        "NG15": ng15_names(),
        "NG15_TIMING": par_tim_rows("NG15_TIMING", loaders.ng15_timing_pairs()),
        "EPTA": par_tim_rows("EPTA", loaders.epta_pairs(epta_variant)),
        "PPTA": par_tim_rows("PPTA", loaders.ppta_pairs()),
        "InPTA_DR1_NB": par_tim_rows("InPTA_DR1_NB", loaders.inpta_pairs("NB")),
        "InPTA_DR1_WB": par_tim_rows("InPTA_DR1_WB", loaders.inpta_pairs("WB")),
        "MPTA_4p5yr": par_tim_rows("MPTA_4p5yr", loaders.mpta_pairs()),
    }


def make_priority_manifest(inventory: dict[str, list[dict]], order: list[str]) -> tuple[list[dict], list[dict]]:
    retained = []
    skipped = []
    used = set()
    for array in order:
        for item in inventory[array]:
            base = item["base_name"]
            if base in used:
                skipped.append(
                    {
                        **item,
                        "status": "skipped",
                        "reason": "duplicate pulsar name skipped by priority-dedup development policy",
                    }
                )
                continue
            used.add(base)
            retained.append({**item, "status": "retained", "reason": "retained by priority-dedup development policy"})
    return retained, skipped


def make_unique_only_manifest(inventory: dict[str, list[dict]], arrays: list[str]) -> tuple[list[dict], list[dict]]:
    counts = Counter(item["base_name"] for array in arrays for item in inventory[array])
    retained = []
    skipped = []
    for array in arrays:
        for item in inventory[array]:
            if counts[item["base_name"]] == 1:
                retained.append({**item, "status": "retained", "reason": "unique across selected arrays"})
            else:
                skipped.append({**item, "status": "skipped", "reason": "duplicate pulsar name excluded by unique-only policy"})
    return retained, skipped


def write_manifest(path: Path, retained: list[dict], skipped: list[dict], policy_name: str, priority_order: list[str]) -> None:
    fields = [
        "index",
        "array",
        "name",
        "base_name",
        "status",
        "reason",
        "source",
        "parfile",
        "timfile",
        "policy",
        "priority_order",
    ]
    rows = []
    for i, item in enumerate(retained):
        rows.append({**item, "index": i, "policy": policy_name, "priority_order": ">".join(priority_order)})
    for item in skipped:
        rows.append({**item, "index": "", "policy": policy_name, "priority_order": ">".join(priority_order)})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def duplicate_table(inventory: dict[str, list[dict]], arrays: list[str]) -> list[dict]:
    by_name: dict[str, list[str]] = defaultdict(list)
    for array in arrays:
        for item in inventory[array]:
            by_name[item["base_name"]].append(array)
    return [
        {"base_name": name, "arrays": sorted(set(arrs)), "occurrences": len(arrs)}
        for name, arrs in sorted(by_name.items())
        if len(set(arrs)) > 1
    ]


def local_count_audit(inventory: dict[str, list[dict]]) -> list[dict]:
    rows: list[dict] = []
    feather_dir = ROOT / "data" / "NG15yr" / "tutorials" / "data" / "feathers"
    feather_names = sorted({path.stem.split("bipm2019-")[-1] for path in feather_dir.glob("*.feather")})
    ng15_release_root = (
        ROOT
        / "参考研究课题-PTA引力波"
        / "data"
        / "nanograv_15yr"
        / "NANOGrav15yr_PulsarTiming_v2.1.0"
    )
    ng15_residual_dir = ng15_release_root / "residuals"
    ng15_reference_names = []
    if ng15_residual_dir.exists():
        ng15_reference_names = sorted(
            {
                path.name.split("_NG15yr_nb.")[0]
                for path in ng15_residual_dir.glob("*_NG15yr_nb.*.res")
            }
        )
    missing_from_feathers = sorted(set(ng15_reference_names) - set(feather_names))
    rows.append(
        {
            "product": "NG15_feathers",
            "status": "LOCAL_COUNT_67_REFERENCE_68",
            "detail": f"Local feather directory contains {len(list(feather_dir.glob('*.feather')))} files; ENTERPRISE loader returns {len(inventory['NG15'])} pulsars.",
            "missing_or_incomplete": ",".join(missing_from_feathers) or "reference inclusion list still required to identify the 68th target",
        }
    )
    if ng15_reference_names:
        j0614_par = ng15_release_root / "narrowband" / "par" / "J0614-3329_PINT_20220301.nb.par"
        j0614_tim = ng15_release_root / "narrowband" / "tim" / "J0614-3329_PINT_20220301.nb.tim"
        ng15_timing_smoke = OUT_DIR / "ng15_timing_loader_smoke.json"
        smoke_detail = "NG15_TIMING loader smoke not yet run."
        if ng15_timing_smoke.exists():
            smoke = json.loads(ng15_timing_smoke.read_text())
            smoke_detail = (
                f"NG15_TIMING loader smoke status={smoke.get('status')}, "
                f"available={smoke.get('available_pulsars')}, "
                f"loaded={len(smoke.get('loaded_pulsars', []))}, "
                f"failures={len(smoke.get('failures', []))}."
            )
        rows.append(
            {
                "product": "NG15_full_timing_reference_tree",
                "status": "REFERENCE_TREE_HAS_68_NAMES_LOADER_SMOKE_PASS",
                "detail": (
                    f"Reference-project NANOGrav 15-year timing residual tree has {len(ng15_reference_names)} unique pulsars. "
                    f"The feather-missing target J0614-3329 has par={j0614_par.exists()} and tim={j0614_tim.exists()}. "
                    f"{smoke_detail}"
                ),
                "missing_or_incomplete": "official NG15 noise-policy wiring remains required before timing-level evidence",
            }
        )

    ppta_meta = ROOT / "data" / "PPTA_DR3" / "_metadata" / "csiro_59381_ppta_dr3__toas_and_parameters__all.json"
    if ppta_meta.exists():
        data = json.loads(ppta_meta.read_text())
        files = [item["filename"] for item in data.get("file", [])]
        tim = sorted({Path(name).stem for name in files if name.endswith(".tim")})
        standard_par = sorted(
            {Path(name).stem for name in files if name.endswith(".par") and not name.endswith("_singlePsrNoise_fit.par")}
        )
        noise_par = sorted({Path(name).stem.replace("_singlePsrNoise_fit", "") for name in files if name.endswith("_singlePsrNoise_fit.par")})
        missing_tim = sorted(set(standard_par) - set(tim))
        rows.append(
            {
                "product": "PPTA_DR3_toas_and_parameters_all",
                "status": "METADATA_HAS_32_PAR_31_TIM",
                "detail": f"Metadata lists {len(standard_par)} standard .par files, {len(noise_par)} single-pulsar-noise .par files, and {len(tim)} .tim files.",
                "missing_or_incomplete": ",".join(missing_tim),
            }
        )
    ppta_github = ROOT / "data" / "PPTA_DR3_github" / "analysis_codes" / "data" / "all"
    if ppta_github.exists():
        tim = sorted({path.stem for path in ppta_github.glob("*.tim")})
        standard_par = sorted({path.stem for path in ppta_github.glob("*.par") if not path.stem.endswith("_singlePsrNoise_fit")})
        head = ""
        git_dir = ROOT / "data" / "PPTA_DR3_github" / ".git"
        if git_dir.exists():
            import subprocess

            head = subprocess.check_output(["git", "-C", str(ROOT / "data" / "PPTA_DR3_github"), "rev-parse", "HEAD"], text=True).strip()
        rows.append(
            {
                "product": "PPTA_DR3_github_analysis_codes_data_all",
                "status": "GITHUB_SOURCE_HAS_32_PAR_32_TIM",
                "detail": f"GitHub source has {len(standard_par)} standard .par files and {len(tim)} .tim files at commit {head}.",
                "missing_or_incomplete": ",".join(sorted(set(standard_par) ^ set(tim))),
            }
        )

    mpta_products = [
        (
            "MPTA_4p5yr_legacy_DOI_partim",
            ROOT / "data" / "MPTA_4p5yr" / "partim" / "noise_subtracted_residuals" / "DOI_partim",
            ROOT / "data" / "MPTA_4p5yr" / "archives" / "partim.tar.gz",
            "legacy local archive previously used by this project",
        ),
        (
            "MPTA_4p5yr_DataCentral_partim_20260425",
            ROOT / "data" / "MPTA_4p5yr" / "partim_datacentral_20260425" / "partim",
            ROOT / "data" / "MPTA_4p5yr" / "archives_redownload_20260425" / "partim_datacentral_20260425.tar.gz",
            "refreshed DataCentral partim download",
        ),
    ]
    for product, mpta_dir, mpta_archive, note in mpta_products:
        if not mpta_dir.exists() and not mpta_archive.exists():
            continue
        local_tim = sorted({path.stem.replace("_16ch", "") for path in mpta_dir.glob("*.tim")})
        local_par = sorted(
            {
                path.stem.removeprefix("MTMSP-").rstrip("-").replace(".par", "")
                for path in mpta_dir.glob("*.par")
            }
        )
        archive_tim_count = ""
        archive_par_count = ""
        if mpta_archive.exists():
            with tarfile.open(mpta_archive, "r:gz") as tar:
                names = tar.getnames()
            archive_tim_count = sum(1 for name in names if name.endswith(".tim"))
            archive_par_count = sum(1 for name in names if name.endswith(".par"))
        pair_count = len(set(local_par) & set(local_tim))
        status = "LOCAL_ARCHIVE_AND_EXTRACTED_TREE_HAVE_83_PAIRS" if pair_count == 83 else "LOCAL_ARCHIVE_AND_EXTRACTED_TREE_HAVE_74_PAIRS"
        rows.append(
            {
                "product": product,
                "status": status,
                "detail": (
                    f"{note}: extracted tree has {len(local_par)} .par names and {len(local_tim)} .tim names "
                    f"({pair_count} matched pairs); archive has {archive_par_count} .par and {archive_tim_count} .tim entries."
                ),
                "missing_or_incomplete": ",".join(sorted(set(local_par) ^ set(local_tim))),
            }
        )

    mpta_supplement = ROOT / "data" / "MPTA_4p5yr" / "archives" / "MPTA_Anisotropy_supplement.zip"
    if mpta_supplement.exists():
        with zipfile.ZipFile(mpta_supplement) as archive:
            names = archive.namelist()
        timing_like = [
            name
            for name in names
            if Path(name).suffix.lower() in {".par", ".tim", ".json", ".csv", ".txt", ".yaml", ".yml"}
            or any(token in Path(name).name.lower() for token in ["manifest", "pulsar", "noise"])
        ]
        mp4_count = sum(1 for name in names if Path(name).suffix.lower() == ".mp4")
        rows.append(
            {
                "product": "MPTA_4p5yr_Anisotropy_supplement_zip",
                "status": "NO_TIMING_INPUTS_FOUND",
                "detail": (
                    f"Supplement zip contains {len(names)} entries, including {mp4_count} MP4 files, "
                    f"and {len(timing_like)} timing-like or manifest-like entries."
                ),
                "missing_or_incomplete": "does not provide timing inputs or duplicate/noise-policy metadata",
            }
        )
    mpta_raw_archive = ROOT / "data" / "MPTA_4p5yr" / "archives" / "archives.tar.gz"
    if mpta_raw_archive.exists():
        with tarfile.open(mpta_raw_archive, "r:gz") as archive:
            members = archive.getmembers()
        names = [member.name for member in members]
        par_count = sum(1 for name in names if Path(name).suffix.lower() == ".par")
        tim_count = sum(1 for name in names if Path(name).suffix.lower() == ".tim")
        dly_count = sum(1 for name in names if Path(name).suffix.lower() == ".dly")
        dirs_count = sum(1 for member in members if member.isdir())
        pulsar_dirs = sorted(
            {
                Path(name).parts[1]
                for name in names
                if len(Path(name).parts) >= 2
                and Path(name).parts[0] == "data_august23_32ch"
                and Path(name).suffix.lower() == ".dly"
            }
        )
        rows.append(
            {
                "product": "MPTA_4p5yr_archives_tar",
                "status": "RAW_DLY_ARCHIVE_NO_PAR_TIM",
                "detail": (
                    f"archives.tar.gz contains {len(members)} entries, {dirs_count} directories, "
                    f"{len(pulsar_dirs)} pulsar directories under data_august23_32ch, {dly_count} .dly files, "
                    f"and {par_count} .par / {tim_count} .tim entries."
                ),
                "missing_or_incomplete": "raw .dly archive does not provide baseline-ready .par/.tim inputs",
            }
        )
    return rows


def gate_items(inventory: dict[str, list[dict]]) -> list[dict]:
    inpta_dr2_root = ROOT / "data" / "InPTA_DR2"
    ref = REFERENCE_BASELINE["dataset_counts"]
    local_ng_feather = len(inventory["NG15"])
    local_ng_timing = len(inventory["NG15_TIMING"])
    local_epta = len(inventory["EPTA"])
    local_ppta = len(inventory["PPTA"])
    local_inpta_nb = len(inventory["InPTA_DR1_NB"])
    local_inpta_wb = len(inventory["InPTA_DR1_WB"])
    local_mpta = len(inventory["MPTA_4p5yr"])
    return [
        {
            "gate": "NG15 timing products",
            "status": "PRESENT_LOCALLY" if local_ng_timing == ref["NG15"] else "COUNT_MISMATCH",
            "detail": (
                f"NANOGrav 15-year reference par/tim target count {local_ng_timing}; reference target count {ref['NG15']}. "
                f"The separate feather development loader has {local_ng_feather} pulsars and must not be mixed into exact timing-level evidence."
            ),
        },
        {
            "gate": "EPTA DR2 timing products",
            "status": "PRESENT_LOCALLY" if local_epta == ref["EPTA"] else "COUNT_MISMATCH",
            "detail": f"EPTA DR2 par/tim count {local_epta}; reference target count {ref['EPTA']}.",
        },
        {
            "gate": "PPTA DR3 timing products",
            "status": "COUNT_MISMATCH" if local_ppta != ref["PPTA"] else "PRESENT_LOCALLY",
            "detail": (
                f"PPTA DR3 local par/tim pair count {local_ppta}; reference target count {ref['PPTA']}.  "
                + (
                    "Count now matches when using the GitHub reference source."
                    if local_ppta == ref["PPTA"]
                    else "The missing/incomplete row must be resolved before exact baseline reproduction."
                )
            ),
        },
        {
            "gate": "InPTA DR1 exact target inputs",
            "status": "PRESENT_BUT_POLICY_BLOCKED" if local_inpta_nb == ref["InPTA_DR1"] else "COUNT_MISMATCH",
            "detail": (
                f"Local tree contains {local_inpta_nb} NB and {local_inpta_wb} WB InPTA DR1 par/tim pairs; "
                f"reference target count is {ref['InPTA_DR1']}.  Read-only diagnostics show NB+tempo2 "
                "loads 13/14 and fails on J0751+1807_NB, WB+PINT loads 13/14 and fails on J0613-0200_WB, "
                "NB+PINT fails all 14 on TCB/TDB handling, and WB+tempo2 fails all 14 with empty-TOA "
                "loader errors.  Exact baseline reproduction therefore remains blocked until the published "
                "NB/WB, timing-package, and TCB-to-TDB preprocessing policy is identified."
            ),
            "local_path": rel(inpta_dr2_root),
            "local_path_exists": inpta_dr2_root.exists(),
        },
        {
            "gate": "MPTA exact target inputs",
            "status": "COUNT_MISMATCH_AND_POLICY_BLOCKED" if local_mpta != ref["MPTA_4p5yr"] else "PRESENT_BUT_POLICY_BLOCKED",
            "detail": (
                f"Local MPTA par/tim pair count is {local_mpta}; reference target count is {ref['MPTA_4p5yr']}. "
                "Production baseline use also requires exact noise-policy and duplicate-combination harmonization."
            ),
        },
        {
            "gate": "duplicate pulsar handling",
            "status": "BLOCKED_BY_POLICY_HARMONIZATION",
            "detail": (
                "Current reduced 3-array diagnostics use priority deduplication or unique-only "
                "dropping.  A production baseline must reproduce the published direct-combination "
                "policy for duplicated pulsars before family evidence."
            ),
        },
        {
            "gate": "timing-level baseline reproduction",
            "status": "NOT_STARTED_FOR_EXACT_TARGET",
            "detail": (
                "No exact public five-PTA baseline likelihood reproduction has passed locally; "
                "family evidence remains on hold."
            ),
        },
    ]


def write_markdown(path: Path, result: dict) -> None:
    lines = [
        "# Manifest v2 and Exact Timing-Level Baseline Gate",
        "",
        f"Generated: {result['generated']}",
        "",
        "This file separates local ingestion status from public availability.  A `BLOCKED`",
        "status below means not yet present locally or not yet harmonized with this",
        "project's exact reproduction manifest; it does not mean the public input is absent.",
        "",
        "## Public Five-PTA Reference Target",
        "",
        f"- Reference source: `{REFERENCE_BASELINE['source']}`.",
        f"- Direct-combination priority: `{REFERENCE_BASELINE['direct_combination_priority']}`.",
        f"- White-noise policy summary: `{REFERENCE_BASELINE['noise_policy']}`.",
        "",
        "| array | reference count | local count |",
        "|---|---:|---:|",
    ]
    local_counts = {row["array"]: row["count"] for row in result["inventory_summary"]}
    reference_rows = [
        ("EPTA", "EPTA"),
        ("InPTA_DR1", "InPTA_DR1_NB"),
        ("MPTA_4p5yr", "MPTA_4p5yr"),
        ("NG15", "NG15_TIMING"),
        ("PPTA", "PPTA"),
    ]
    for ref_name, local_name in reference_rows:
        lines.append(f"| `{ref_name}` | {REFERENCE_BASELINE['dataset_counts'][ref_name]} | {local_counts.get(local_name, 0)} |")

    lines.extend(
        [
            "",
            "## Local Inventory",
            "",
            "| product | count | status |",
            "|---|---:|---|",
        ]
    )
    for row in result["inventory_summary"]:
        lines.append(f"| `{row['array']}` | {row['count']} | `{row['status']}` |")

    lines.extend(
        [
            "",
            "## Reduced Development Namespace",
            "",
            "The current three-array timing diagnostics use an `NG15-priority deduplicated",
            "reduced 3-array development namespace`, not a full public five-PTA direct",
            "combination.",
            "",
            f"- Retained pulsars: `{result['ng15_priority_manifest']['retained']}`.",
            f"- Skipped duplicate occurrences: `{result['ng15_priority_manifest']['skipped']}`.",
            f"- Manifest CSV: `{result['ng15_priority_manifest']['csv']}`.",
            f"- Unique-only retained pulsars: `{result['unique_only_manifest']['retained']}`.",
            f"- Unique-only skipped duplicate occurrences: `{result['unique_only_manifest']['skipped']}`.",
            f"- Unique-only CSV: `{result['unique_only_manifest']['csv']}`.",
            "",
            "## Duplicate Sources Across NG15/EPTA/PPTA",
            "",
            "| pulsar | arrays | occurrences |",
            "|---|---|---:|",
        ]
    )
    for row in result["duplicates_ng15_epta_ppta"]:
        lines.append(f"| `{row['base_name']}` | `{', '.join(row['arrays'])}` | {row['occurrences']} |")

    lines.extend(
        [
            "",
            "## Local Count Discrepancy Audit",
            "",
            "| product | status | detail | missing/incomplete |",
            "|---|---|---|---|",
        ]
    )
    for row in result["local_count_audit"]:
        lines.append(
            f"| `{row['product']}` | `{row['status']}` | {row['detail']} | `{row.get('missing_or_incomplete', '')}` |"
        )

    lines.extend(["", "## Exact Baseline Gate", "", "| gate | status | detail |", "|---|---|---|"])
    for row in result["gate_items"]:
        lines.append(f"| {row['gate']} | `{row['status']}` | {row['detail']} |")

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Do not extend the current reduced free-red chains as a production route.",
            "- First reproduce the selected public timing-level baseline with the same duplicate-handling and noise policy.",
            "- Only after that gate passes should minimal family evidence be run on timing-level likelihoods.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epta-variant", default="DR2new+", choices=["DR2new+", "DR2new", "DR2full+", "DR2full"])
    parser.add_argument("--tag", default="2026-04-24")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    THEORY_DIR.mkdir(parents=True, exist_ok=True)

    inv = local_inventory(args.epta_variant)
    selected_arrays = ["NG15", "EPTA", "PPTA"]
    retained, skipped = make_priority_manifest(inv, selected_arrays)
    unique_retained, unique_skipped = make_unique_only_manifest(inv, selected_arrays)

    manifest_csv = OUT_DIR / f"manifest_v2_ng15_priority_dedup_3array_development_{args.tag}.csv"
    unique_csv = OUT_DIR / f"manifest_v2_unique_only_3array_development_{args.tag}.csv"
    write_manifest(manifest_csv, retained, skipped, "ng15_priority_deduplicated_reduced_3array_development", selected_arrays)
    write_manifest(unique_csv, unique_retained, unique_skipped, "unique_only_reduced_3array_development", [])

    result = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "status": "BLOCKED_BEFORE_TIMING_LEVEL_FAMILY_EVIDENCE",
        "epta_variant": args.epta_variant,
        "inventory_summary": [
            {"array": array, "count": len(rows), "status": "PRESENT_LOCALLY" if rows else "NOT_PRESENT_LOCALLY"}
            for array, rows in inv.items()
        ],
        "ng15_priority_manifest": {"retained": len(retained), "skipped": len(skipped), "csv": rel(manifest_csv)},
        "unique_only_manifest": {"retained": len(unique_retained), "skipped": len(unique_skipped), "csv": rel(unique_csv)},
        "duplicates_ng15_epta_ppta": duplicate_table(inv, selected_arrays),
        "local_count_audit": local_count_audit(inv),
        "gate_items": gate_items(inv),
    }

    json_path = OUT_DIR / f"manifest_v2_exact_baseline_gate_{args.tag}.json"
    md_path = THEORY_DIR / f"PTA_manifest_v2_exact_baseline_gate_{args.tag}.md"
    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    write_markdown(md_path, result)
    print(json.dumps({"status": result["status"], "json": rel(json_path), "markdown": rel(md_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
