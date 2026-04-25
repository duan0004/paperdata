#!/usr/bin/env python3
"""LSS tomography production gate for 2MPZ, WISExSCOS, and DESI.

This gate deliberately separates template production from science claims:

1. stage a compact 2MPZ RA/Dec/photo-z extract from the official SSA SQL
   endpoint, because the previous all-sky direct-cone attempt is an error file;
2. build redshift-binned 2MPZ and WISExSCOS HEALPix overdensity maps;
3. project 2MPZ and WISExSCOS maps onto the NG15 pulsar-pair ORF geometry;
4. inventory the already-built DESI redshift-binned ORFs;
5. mark the NG15+2MPZ published-null reproduction as blocked until a
   selection-corrected 2MPZ likelihood reproduces the reference null.

The output is a production gate, not an anisotropy detection and not a
Bayes-factor calculation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import re
import shutil
import sys
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = ROOT / "code"
DATA = ROOT / "data" / "LSS"
MAP_DIR = DATA / "maps"
ORF_DIR = DATA / "orf_templates"
OUT_DIR = ROOT / "results" / "lss_tomography"

if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import lss_orf_null_tests as lss_orf  # noqa: E402


@dataclass(frozen=True)
class RedshiftBin:
    key: str
    z_min: float | None
    z_max: float | None
    include_upper: bool = False

    def contains(self, z: np.ndarray) -> np.ndarray:
        mask = np.isfinite(z)
        if self.z_min is not None:
            mask &= z >= self.z_min
        if self.z_max is not None:
            if self.include_upper:
                mask &= z <= self.z_max
            else:
                mask &= z < self.z_max
        return mask


BINS = [
    RedshiftBin("all", None, None),
    RedshiftBin("z_0_0p1", 0.0, 0.1),
    RedshiftBin("z_0p1_0p2", 0.1, 0.2),
    RedshiftBin("z_0_0p2", 0.0, 0.2),
    RedshiftBin("z_0p2_0p3", 0.2, 0.3),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nside", type=int, default=64)
    parser.add_argument("--tag", default="prod_nside64_20260424")
    parser.add_argument("--nulls", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260424)
    parser.add_argument("--skip-download", action="store_true")
    return parser.parse_args()


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


def download_2mpz_compact(path: Path, skip_download: bool) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "catalog": "2MPZ",
        "path": rel(path),
        "status": "pending",
        "bytes": "",
        "query": "SELECT ra,dec,zPhoto FROM TWOMPZ..twompzPhotoz",
        "error": "",
    }
    if path.exists() and path.stat().st_size > 0:
        row["status"] = "exists"
        row["bytes"] = path.stat().st_size
        return row
    if skip_download:
        row["status"] = "missing_skip_download"
        return row

    endpoint = "http://ssa.roe.ac.uk:8080/ssa/SSASQL"
    payload = urllib.parse.urlencode(
        {
            "action": "freeform",
            "server": "amachine",
            "sqlstmt": row["query"],
            "format": "CSV",
            "compress": "GZIP",
            "rows": "30",
            "emailAddress": "",
        }
    ).encode()
    tmp = path.with_suffix(path.suffix + ".part")
    try:
        req = urllib.request.Request(endpoint, data=payload, headers={"User-Agent": "Codex-PTA-LSS-gate/1.0"})
        html = urllib.request.urlopen(req, timeout=180).read().decode("ISO-8859-1", "ignore")
        match = re.search(r'href="(http://ssa\.roe\.ac\.uk/tmp/[^"]+\.csv\.gz)"', html)
        if not match:
            raise RuntimeError("SSA SQL response did not contain a downloadable .csv.gz result link")
        with urllib.request.urlopen(match.group(1), timeout=300) as src, tmp.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        tmp.replace(path)
        row["status"] = "downloaded"
        row["bytes"] = path.stat().st_size
        row["download_url"] = match.group(1)
    except Exception as exc:
        if tmp.exists():
            tmp.unlink()
        row["status"] = "failed"
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def make_empty_maps(nside: int) -> dict[str, np.ndarray]:
    import healpy as hp

    return {zbin.key: np.zeros(hp.nside2npix(nside), dtype=np.float64) for zbin in BINS}


def build_2mpz_maps(path: Path, nside: int, tag: str) -> tuple[list[dict], list[dict]]:
    import healpy as hp

    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    maps = make_empty_maps(nside)
    rows_read = 0
    rows_valid = 0
    z_min = math.inf
    z_max = -math.inf
    with gzip.open(path, "rt", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows_read += 1
            try:
                ra = float(row["ra"])
                dec = float(row["dec"])
                z = float(row["zPhoto"])
            except Exception:
                continue
            if not (math.isfinite(ra) and math.isfinite(dec) and math.isfinite(z) and z >= 0.0):
                continue
            pix = hp.ang2pix(nside, ra, dec, lonlat=True)
            rows_valid += 1
            z_min = min(z_min, z)
            z_max = max(z_max, z)
            for zbin in BINS:
                if zbin.contains(np.asarray([z]))[0]:
                    maps[zbin.key][pix] += 1.0
    summary_rows = save_overdensity_maps("twompz", maps, nside, tag, selection_mask=None)
    input_rows = [
        {
            "catalog": "2MPZ",
            "path": rel(path),
            "rows_read": rows_read,
            "rows_valid": rows_valid,
            "z_min": z_min if math.isfinite(z_min) else "",
            "z_max": z_max if math.isfinite(z_max) else "",
            "selection_note": "No random catalog applied; this is a template-production gate only.",
        }
    ]
    return summary_rows, input_rows


def degrade_wise_mask(nside: int) -> np.ndarray:
    import healpy as hp

    mask_path = DATA / "WISExSCOS" / "WISExSCOSmask.fits.gz"
    if not mask_path.exists():
        raise FileNotFoundError(mask_path)
    mask = np.asarray(hp.read_map(str(mask_path), verbose=False), dtype=np.float64)
    if hp.npix2nside(mask.size) != nside:
        mask = hp.ud_grade(mask, nside, order_in="RING", order_out="RING", power=0)
    return mask > 0.5


def build_wisexscos_maps(path: Path, nside: int, tag: str) -> tuple[list[dict], list[dict]]:
    import healpy as hp

    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    maps = make_empty_maps(nside)
    selection_mask = degrade_wise_mask(nside)
    rows_read = 0
    rows_valid = 0
    z_min = math.inf
    z_max = -math.inf
    with gzip.open(path, "rt", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows_read += 1
            try:
                ra = float(row["ra"])
                dec = float(row["dec"])
                z = float(row["zPhoto_Corr"] or row["zPhoto_ANN"])
            except Exception:
                continue
            if not (math.isfinite(ra) and math.isfinite(dec) and math.isfinite(z) and z >= 0.0):
                continue
            pix = hp.ang2pix(nside, ra, dec, lonlat=True)
            if not selection_mask[pix]:
                continue
            rows_valid += 1
            z_min = min(z_min, z)
            z_max = max(z_max, z)
            for zbin in BINS:
                if zbin.contains(np.asarray([z]))[0]:
                    maps[zbin.key][pix] += 1.0
    summary_rows = save_overdensity_maps("wisexscos", maps, nside, tag, selection_mask=selection_mask)
    input_rows = [
        {
            "catalog": "WISExSCOS",
            "path": rel(path),
            "rows_read": rows_read,
            "rows_valid_after_mask": rows_valid,
            "z_min": z_min if math.isfinite(z_min) else "",
            "z_max": z_max if math.isfinite(z_max) else "",
            "selection_note": "Uses the public binary WISExSCOS mask degraded to the production nside.",
        }
    ]
    return summary_rows, input_rows


def save_overdensity_maps(
    prefix: str,
    maps: dict[str, np.ndarray],
    nside: int,
    tag: str,
    selection_mask: np.ndarray | None,
) -> list[dict]:
    MAP_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    npix = next(iter(maps.values())).size
    for zbin in BINS:
        key = zbin.key
        counts = maps[key]
        valid = np.ones(npix, dtype=bool) if selection_mask is None else selection_mask.copy()
        valid &= np.isfinite(counts)
        mean_count = float(np.mean(counts[valid])) if np.any(valid) else math.nan
        delta = np.zeros(npix, dtype=np.float64)
        if math.isfinite(mean_count) and mean_count > 0.0:
            delta[valid] = counts[valid] / mean_count - 1.0
        else:
            valid[:] = False
        stem = MAP_DIR / f"{prefix}_{key}_{tag}"
        counts_path = stem.with_name(stem.name + "_counts.npy")
        delta_path = stem.with_name(stem.name + "_delta.npy")
        mask_path = stem.with_name(stem.name + "_mask.npy")
        np.save(counts_path, counts)
        np.save(delta_path, delta)
        np.save(mask_path, valid.astype(np.uint8))
        rows.append(
            {
                "catalog": prefix,
                "bin": key,
                "nside": nside,
                "counts_sum": float(np.sum(counts)),
                "mean_count_valid": mean_count,
                "valid_pixel_fraction": float(np.mean(valid)),
                "occupied_pixel_fraction_valid": float(np.mean((counts > 0.0) & valid)) if np.any(valid) else "",
                "delta_std_valid": float(np.std(delta[valid])) if np.any(valid) else "",
                "counts_map": rel(counts_path),
                "delta_map": rel(delta_path),
                "mask_map": rel(mask_path),
            }
        )
    return rows


def prepare_field(delta_path: Path, mask_path: Path, nside: int) -> tuple[np.ndarray, dict]:
    import healpy as hp

    delta = np.asarray(np.load(delta_path), dtype=np.float64)
    mask = np.asarray(np.load(mask_path), dtype=bool)
    npix = hp.nside2npix(nside)
    if delta.size != npix or mask.size != npix:
        raise RuntimeError(f"npix mismatch for {delta_path}")
    valid = mask & np.isfinite(delta)
    if not np.any(valid):
        raise RuntimeError(f"empty valid mask for {delta_path}")
    field = np.zeros(npix, dtype=np.float64)
    field[valid] = delta[valid] - float(np.mean(delta[valid]))
    rms = float(np.sqrt(np.mean(field[valid] ** 2)))
    if rms > 0.0:
        field[valid] /= rms
    return field, {"valid_pixels": int(np.count_nonzero(valid)), "mask_fraction": float(np.mean(valid)), "rms_before_norm": rms}


def build_orfs_for_catalog(catalog: str, map_rows: list[dict], nside: int, tag: str, nulls: int, rng: np.random.Generator) -> tuple[list[dict], list[dict]]:
    import healpy as hp

    ORF_DIR.mkdir(parents=True, exist_ok=True)
    names, pos = lss_orf.load_ng15_positions()
    pairs, hd, pair_rows = lss_orf.pairs_and_hd(names, pos)
    pair_path = ORF_DIR / f"ng15_pair_index_{tag}.csv"
    if not pair_path.exists():
        write_csv(pair_path, pair_rows, ["i", "j", "pulsar_i", "pulsar_j", "cos_zeta"])
    npix = hp.nside2npix(nside)
    omega = np.asarray(hp.pix2vec(nside, np.arange(npix))).T
    fplus, fcross = lss_orf.antenna_response(pos, omega)

    orf_rows = []
    null_rows = []
    for row in map_rows:
        delta_path = ROOT / row["delta_map"]
        mask_path = ROOT / row["mask_map"]
        field, field_summary = prepare_field(delta_path, mask_path, nside)
        vec = lss_orf.pair_vector(field, fplus, fcross, pairs)
        stats = lss_orf.vector_stats(vec, hd)
        out_path = ORF_DIR / f"{catalog}_{row['bin']}_{tag}_ng15_orf_pairs.npy"
        np.save(out_path, vec)
        orf_rows.append(
            {
                "catalog": catalog,
                "bin": row["bin"],
                **field_summary,
                **stats,
                "orf_vector": rel(out_path),
                "pair_index": rel(pair_path),
            }
        )
        null_rows.append({"catalog": catalog, "bin": row["bin"], **random_map_null(field, fplus, fcross, pairs, hd, stats, nulls, rng)})
    return orf_rows, null_rows


def random_map_null(
    field: np.ndarray,
    fplus: np.ndarray,
    fcross: np.ndarray,
    pairs: list[tuple[int, int]],
    hd: np.ndarray,
    real_stats: dict,
    n_nulls: int,
    rng: np.random.Generator,
) -> dict:
    valid = field != 0.0
    values = field[valid].copy()
    hd_corrs = []
    l2s = []
    for _ in range(n_nulls):
        shuffled = np.zeros_like(field)
        shuffled_values = values.copy()
        rng.shuffle(shuffled_values)
        shuffled[valid] = shuffled_values
        vec = lss_orf.pair_vector(shuffled, fplus, fcross, pairs)
        stats = lss_orf.vector_stats(vec, hd)
        hd_corrs.append(stats["pearson_with_hd"])
        l2s.append(stats["pair_l2"])
    hd_corrs = np.asarray(hd_corrs, dtype=np.float64)
    l2s = np.asarray(l2s, dtype=np.float64)
    real_corr = float(real_stats["pearson_with_hd"])
    return {
        "n_nulls": n_nulls,
        "real_hd_corr": real_corr,
        "null_hd_corr_mean": float(np.nanmean(hd_corrs)),
        "null_hd_corr_std": float(np.nanstd(hd_corrs)),
        "real_hd_corr_z": float((real_corr - np.nanmean(hd_corrs)) / max(np.nanstd(hd_corrs), 1.0e-30)),
        "real_abs_hd_corr_percentile": float(np.mean(np.abs(hd_corrs) <= abs(real_corr))),
        "real_l2_z": float((float(real_stats["pair_l2"]) - np.mean(l2s)) / max(np.std(l2s), 1.0e-30)),
    }


def inventory_desi_orfs() -> list[dict]:
    rows = []
    for path in sorted(ORF_DIR.glob("desi_bgs_bright21p5_*_prod_nside64_weight_nrand36_ng15_orf_pairs.npy")):
        arr = np.asarray(np.load(path), dtype=np.float64)
        rows.append(
            {
                "catalog": "DESI_BGS_BRIGHT21p5",
                "orf_vector": rel(path),
                "pair_count": int(arr.size),
                "pair_mean": float(np.mean(arr)),
                "pair_std": float(np.std(arr)),
                "status": "present",
            }
        )
    return rows


def write_markdown(result: dict, download_rows: list[dict], map_rows: list[dict], orf_rows: list[dict], null_rows: list[dict], desi_rows: list[dict]) -> None:
    lines = [
        "# LSS Tomography Production Gate",
        "",
        f"Generated: {result['generated']}",
        "",
        "## Decision Boundary",
        "",
        "- NG15+2MPZ published-null reproduction is the first science gate.",
        "- DESI/WISExSCOS redshift-binned ORFs are production inputs only until that null is reproduced.",
        "- No positive LSS-correlated signal should be trusted or reported from this gate.",
        "",
        f"Overall status: `{result['status']}`.",
        f"2MPZ null reproduction status: `{result['ng15_2mpz_null_reproduction']}`.",
        "",
        "## 2MPZ Download",
        "",
        "| catalog | status | bytes | path | query/error |",
        "|---|---|---:|---|---|",
    ]
    for row in download_rows:
        msg = row.get("error") or row.get("query", "")
        lines.append(f"| {row['catalog']} | `{row['status']}` | {row.get('bytes', '')} | `{row['path']}` | {msg} |")

    lines.extend(["", "## Map Products", "", "| catalog | bin | counts | valid fraction | delta std | delta map |", "|---|---|---:|---:|---:|---|"])
    for row in map_rows:
        lines.append(
            f"| `{row['catalog']}` | `{row['bin']}` | {float(row['counts_sum']):.6g} | "
            f"{float(row['valid_pixel_fraction']):.6f} | {float(row['delta_std_valid']) if row['delta_std_valid'] != '' else float('nan'):.6g} | `{row['delta_map']}` |"
        )

    lines.extend(["", "## NG15 ORF Templates", "", "| catalog | bin | corr with HD | null z | abs-percentile | ORF vector |", "|---|---|---:|---:|---:|---|"])
    null_by_key = {(row["catalog"], row["bin"]): row for row in null_rows}
    for row in orf_rows:
        null = null_by_key[(row["catalog"], row["bin"])]
        lines.append(
            f"| `{row['catalog']}` | `{row['bin']}` | {float(row['pearson_with_hd']):.6f} | "
            f"{float(null['real_hd_corr_z']):.3f} | {float(null['real_abs_hd_corr_percentile']):.3f} | `{row['orf_vector']}` |"
        )

    lines.extend(["", "## DESI ORF Inventory", "", "| status | pair count | pair std | ORF vector |", "|---|---:|---:|---|"])
    for row in desi_rows:
        lines.append(f"| `{row['status']}` | {row['pair_count']} | {float(row['pair_std']):.6g} | `{row['orf_vector']}` |")

    lines.extend(
        [
            "",
            "## Blockers Before Science Use",
            "",
            "1. Reproduce the published NG15+2MPZ null with the same 2MPZ template, selection/mask treatment, and likelihood definition.",
            "2. Add isotropic and random-map null injections inside the ENTERPRISE anisotropy likelihood.",
            "3. Run posterior sampling/evidence only after the null and injection gates pass.",
            "",
        ]
    )
    (OUT_DIR / f"lss_tomography_production_gate_{result['tag']}.md").write_text("\n".join(lines))


def main() -> int:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MAP_DIR.mkdir(parents=True, exist_ok=True)
    ORF_DIR.mkdir(parents=True, exist_ok=True)

    download_rows = []
    map_rows = []
    input_rows = []
    orf_rows = []
    null_rows = []

    twompz_path = DATA / "2MPZ" / "twompz_ra_dec_zphoto.csv.gz"
    row = download_2mpz_compact(twompz_path, args.skip_download)
    download_rows.append(row)

    if row["status"] in {"exists", "downloaded"}:
        rows, inputs = build_2mpz_maps(twompz_path, args.nside, args.tag)
        map_rows.extend(rows)
        input_rows.extend(inputs)
    else:
        raise RuntimeError(f"2MPZ compact download is required before map/ORF production: {row}")

    wise_path = DATA / "WISExSCOS" / "wiseScosPhotoz160708.csv.gz"
    wise_rows, wise_inputs = build_wisexscos_maps(wise_path, args.nside, args.tag)
    map_rows.extend(wise_rows)
    input_rows.extend(wise_inputs)

    rng = np.random.default_rng(args.seed)
    for catalog in ("twompz", "wisexscos"):
        selected = [row for row in map_rows if row["catalog"] == catalog]
        rows, nulls = build_orfs_for_catalog(catalog, selected, args.nside, args.tag, args.nulls, rng)
        orf_rows.extend(rows)
        null_rows.extend(nulls)

    desi_rows = inventory_desi_orfs()
    result = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "tag": args.tag,
        "status": "BLOCKED_BEFORE_LSS_SCIENCE_CLAIM",
        "ng15_2mpz_null_reproduction": "NOT_REPRODUCED_YET",
        "nside": args.nside,
        "nulls": args.nulls,
        "download_manifest": rel(OUT_DIR / f"lss_tomography_production_downloads_{args.tag}.csv"),
        "map_summary": rel(OUT_DIR / f"lss_tomography_production_maps_{args.tag}.csv"),
        "input_summary": rel(OUT_DIR / f"lss_tomography_production_inputs_{args.tag}.csv"),
        "orf_summary": rel(OUT_DIR / f"lss_tomography_production_orfs_{args.tag}.csv"),
        "null_summary": rel(OUT_DIR / f"lss_tomography_production_nulls_{args.tag}.csv"),
        "desi_inventory": rel(OUT_DIR / f"lss_tomography_production_desi_inventory_{args.tag}.csv"),
        "report": rel(OUT_DIR / f"lss_tomography_production_gate_{args.tag}.md"),
    }

    write_csv(OUT_DIR / f"lss_tomography_production_downloads_{args.tag}.csv", download_rows, ["catalog", "status", "bytes", "path", "query", "download_url", "error"])
    write_csv(OUT_DIR / f"lss_tomography_production_maps_{args.tag}.csv", map_rows, ["catalog", "bin", "nside", "counts_sum", "mean_count_valid", "valid_pixel_fraction", "occupied_pixel_fraction_valid", "delta_std_valid", "counts_map", "delta_map", "mask_map"])
    write_csv(OUT_DIR / f"lss_tomography_production_inputs_{args.tag}.csv", input_rows, sorted({k for row in input_rows for k in row}))
    write_csv(OUT_DIR / f"lss_tomography_production_orfs_{args.tag}.csv", orf_rows, ["catalog", "bin", "valid_pixels", "mask_fraction", "rms_before_norm", "pair_mean", "pair_std", "pair_l2", "pearson_with_hd", "max_abs_pair", "orf_vector", "pair_index"])
    write_csv(OUT_DIR / f"lss_tomography_production_nulls_{args.tag}.csv", null_rows, ["catalog", "bin", "n_nulls", "real_hd_corr", "null_hd_corr_mean", "null_hd_corr_std", "real_hd_corr_z", "real_abs_hd_corr_percentile", "real_l2_z"])
    write_csv(OUT_DIR / f"lss_tomography_production_desi_inventory_{args.tag}.csv", desi_rows, ["catalog", "status", "pair_count", "pair_mean", "pair_std", "orf_vector"])
    (OUT_DIR / f"lss_tomography_production_gate_{args.tag}.json").write_text(
        json.dumps(
            {
                "result": result,
                "downloads": download_rows,
                "inputs": input_rows,
                "maps": map_rows,
                "orfs": orf_rows,
                "nulls": null_rows,
                "desi_inventory": desi_rows,
            },
            indent=2,
            ensure_ascii=False,
            sort_keys=True,
        )
        + "\n"
    )
    write_markdown(result, download_rows, map_rows, orf_rows, null_rows, desi_rows)
    print(json.dumps({"status": result["status"], "report": result["report"], "json": rel(OUT_DIR / f"lss_tomography_production_gate_{args.tag}.json")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
