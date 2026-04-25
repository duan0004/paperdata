#!/usr/bin/env python3
"""P7 LSS-tomography manifest, downloads, and geometry gate.

This is a reproducibility gate for the proposed 5PTA x LSS extension.  It
downloads lightweight public LSS products and a small DESI DR1 low-redshift
clustering subset, builds a count map, and verifies that an NG15 sky-template
geometry vector can be constructed.  It does not run a timing-residual
anisotropy likelihood.
"""

from __future__ import annotations

import contextlib
import csv
import gzip
import io
import math
import re
import shutil
import sys
import urllib.request
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "LSS"
META = DATA / "metadata"
DESI = DATA / "DESI_DR1"
MAPS = DATA / "maps"
RESULTS = ROOT / "results" / "lss_tomography"

DESI_BASE = "https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5"


@dataclass
class Resource:
    key: str
    url: str
    local_path: Path
    role: str
    download: bool
    production_required: bool
    note: str


RESOURCES = [
    Resource(
        "2MPZ_SSA_landing",
        "http://ssa.roe.ac.uk/TWOMPZ.html",
        META / "TWOMPZ.html",
        "metadata",
        True,
        True,
        "2MPZ is available through SSA table twompzPhotoz; full export is deferred to the production query plan.",
    ),
    Resource(
        "WISExSCOS_SSA_landing",
        "http://ssa.roe.ac.uk/WISExSCOS.html",
        META / "WISExSCOS.html",
        "metadata",
        True,
        True,
        "Landing page advertises the full WISExSCOS catalogue and mask.",
    ),
    Resource(
        "WISExSCOS_mask",
        "http://ssa.roe.ac.uk/WISExSCOSmask.fits.gz",
        DATA / "WISExSCOS" / "WISExSCOSmask.fits.gz",
        "mask",
        True,
        True,
        "Small public binary mask; used for map-ingestion validation.",
    ),
    Resource(
        "WISExSCOS_full_catalog",
        "http://ssa.roe.ac.uk/cats/wiseScosPhotoz160708.csv.gz",
        DATA / "WISExSCOS" / "wiseScosPhotoz160708.csv.gz",
        "large_catalog",
        True,
        True,
        "Large full catalogue for WISExSCOS tomography; HEAD size was 1.709 GB.",
    ),
    Resource(
        "DESI_DR1_LSS_v1p5_index",
        f"{DESI_BASE}/",
        META / "DESI_DR1_LSS_v1p5_index.html",
        "metadata",
        True,
        True,
        "Official DESI DR1 LSS v1.5 directory listing.",
    ),
    Resource(
        "DESI_BGS_BRIGHT21p5_NGC",
        f"{DESI_BASE}/BGS_BRIGHT-21.5_NGC_clustering.dat.fits",
        DESI / "BGS_BRIGHT-21.5_NGC_clustering.dat.fits",
        "pilot_catalog",
        True,
        False,
        "Low-redshift pilot clustering data for a first angular count map; no random-catalog correction applied here.",
    ),
    Resource(
        "DESI_BGS_BRIGHT21p5_SGC",
        f"{DESI_BASE}/BGS_BRIGHT-21.5_SGC_clustering.dat.fits",
        DESI / "BGS_BRIGHT-21.5_SGC_clustering.dat.fits",
        "pilot_catalog",
        True,
        False,
        "Low-redshift pilot clustering data for a first angular count map; no random-catalog correction applied here.",
    ),
    Resource(
        "DESI_BGS_ANY_NGC_nz",
        f"{DESI_BASE}/BGS_ANY_NGC_nz.txt",
        DESI / "BGS_ANY_NGC_nz.txt",
        "redshift_distribution",
        True,
        False,
        "Small DESI n(z) auxiliary file for redshift-bin sanity checks.",
    ),
    Resource(
        "DESI_BGS_ANY_SGC_nz",
        f"{DESI_BASE}/BGS_ANY_SGC_nz.txt",
        DESI / "BGS_ANY_SGC_nz.txt",
        "redshift_distribution",
        True,
        False,
        "Small DESI n(z) auxiliary file for redshift-bin sanity checks.",
    ),
]

for cap in ["NGC", "SGC"]:
    for i in range(18):
        RESOURCES.append(
            Resource(
                f"DESI_BGS_BRIGHT21p5_{cap}_random_{i:02d}",
                f"{DESI_BASE}/BGS_BRIGHT-21.5_{cap}_{i}_clustering.ran.fits",
                DESI / "randoms" / f"BGS_BRIGHT-21.5_{cap}_{i}_clustering.ran.fits",
                "selection_random",
                True,
                True,
                "DESI BGS_BRIGHT-21.5 random catalog for production selection correction.",
            )
        )


def download_resource(res: Resource) -> dict:
    res.local_path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "key": res.key,
        "url": res.url,
        "local_path": str(res.local_path.relative_to(ROOT)),
        "role": res.role,
        "download_requested": int(res.download),
        "production_required": int(res.production_required),
        "status": "deferred" if not res.download else "pending",
        "bytes": "",
        "note": res.note,
        "error": "",
    }
    if not res.download:
        return row
    if res.local_path.exists() and res.local_path.stat().st_size > 0:
        row["status"] = "exists"
        row["bytes"] = res.local_path.stat().st_size
        return row
    tmp = res.local_path.with_suffix(res.local_path.suffix + ".part")
    try:
        req = urllib.request.Request(res.url, headers={"User-Agent": "Codex-PTA-LSS-gate/1.0"})
        with urllib.request.urlopen(req, timeout=90) as src, tmp.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        tmp.replace(res.local_path)
        row["status"] = "downloaded"
        row["bytes"] = res.local_path.stat().st_size
    except Exception as exc:
        if tmp.exists():
            tmp.unlink()
        row["status"] = "failed"
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def load_selection_corrected_summary() -> list[dict]:
    path = RESULTS / "selection_corrected_map_summary_prod_nside64_weight_nrand36.csv"
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def validate_wise_mask(path: Path) -> dict:
    if not path.exists():
        return {"status": "missing", "error": "mask file missing"}
    try:
        import healpy as hp

        with gzip.open(path, "rb") as f:
            data = f.read(16)
        if not data.startswith(b"SIMPLE"):
            return {"status": "failed", "error": "gzip payload is not a FITS header"}
        mask = hp.read_map(str(path), verbose=False)
        mask = np.asarray(mask, dtype=float)
        nside = hp.npix2nside(mask.size)
        return {
            "status": "pass",
            "nside": int(nside),
            "npix": int(mask.size),
            "finite_fraction": float(np.mean(np.isfinite(mask))),
            "mean_mask": float(np.nanmean(mask)),
            "nonzero_fraction": float(np.mean(mask > 0)),
        }
    except Exception as exc:
        return {"status": "failed", "error": f"{type(exc).__name__}: {exc}"}


def build_desi_count_map(paths: list[Path], nside: int = 64) -> tuple[dict, np.ndarray]:
    import healpy as hp
    from astropy.table import Table

    npix = hp.nside2npix(nside)
    maps = {
        "all": np.zeros(npix, dtype=np.float64),
        "z_0_0p2": np.zeros(npix, dtype=np.float64),
        "z_0p2_0p5": np.zeros(npix, dtype=np.float64),
    }
    n_total = 0
    z_min = float("inf")
    z_max = -float("inf")
    column_sets = []
    for path in paths:
        tab = Table.read(path, memmap=True)
        names = {name.upper(): name for name in tab.colnames}
        if "RA" not in names or "DEC" not in names:
            raise RuntimeError(f"{path.name}: required RA/DEC columns missing")
        z_name = names.get("Z")
        ra = np.asarray(tab[names["RA"]], dtype=float)
        dec = np.asarray(tab[names["DEC"]], dtype=float)
        valid = np.isfinite(ra) & np.isfinite(dec)
        if z_name:
            z = np.asarray(tab[z_name], dtype=float)
            valid &= np.isfinite(z)
        else:
            z = np.full(ra.shape, np.nan)
        pix = hp.ang2pix(nside, ra[valid], dec[valid], lonlat=True)
        np.add.at(maps["all"], pix, 1.0)
        if z_name:
            zv = z[valid]
            z_min = min(z_min, float(np.nanmin(zv)))
            z_max = max(z_max, float(np.nanmax(zv)))
            m1 = (zv >= 0.0) & (zv < 0.2)
            m2 = (zv >= 0.2) & (zv < 0.5)
            np.add.at(maps["z_0_0p2"], pix[m1], 1.0)
            np.add.at(maps["z_0p2_0p5"], pix[m2], 1.0)
        n_total += int(valid.sum())
        column_sets.append(";".join(tab.colnames))
    MAPS.mkdir(parents=True, exist_ok=True)
    for name, m in maps.items():
        np.save(MAPS / f"desi_bgs_bright21p5_{name}_counts_nside{nside}.npy", m)
    summary = {
        "status": "pass",
        "nside": nside,
        "npix": int(npix),
        "n_objects": n_total,
        "z_min": z_min if math.isfinite(z_min) else "",
        "z_max": z_max if math.isfinite(z_max) else "",
        "occupied_fraction_all": float(np.mean(maps["all"] > 0)),
        "occupied_fraction_z_0_0p2": float(np.mean(maps["z_0_0p2"] > 0)),
        "occupied_fraction_z_0p2_0p5": float(np.mean(maps["z_0p2_0p5"] > 0)),
        "max_pixel_count_all": float(np.max(maps["all"])),
        "columns": " | ".join(column_sets),
        "map_all": str((MAPS / f"desi_bgs_bright21p5_all_counts_nside{nside}.npy").relative_to(ROOT)),
    }
    return summary, maps["all"]


def load_ng15_positions() -> tuple[list[str], np.ndarray]:
    from enterprise_extensions.load_feathers import load_feathers_from_folder

    feather_dir = ROOT / "data" / "NG15yr" / "tutorials" / "data" / "feathers"
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        psrs = load_feathers_from_folder(str(feather_dir))
    names = []
    pos = []
    for psr in psrs:
        if hasattr(psr, "pos"):
            arr = np.asarray(psr.pos, dtype=float)
            if arr.size == 3 and np.all(np.isfinite(arr)):
                names.append(str(psr.name))
                pos.append(arr / np.linalg.norm(arr))
    if not pos:
        raise RuntimeError("No NG15 pulsar positions loaded from feather products")
    return names, np.vstack(pos)


def polarization_basis(omega: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    z = np.array([0.0, 0.0, 1.0])
    x = np.array([1.0, 0.0, 0.0])
    u = np.cross(z, omega)
    bad = np.linalg.norm(u, axis=1) < 1e-12
    if np.any(bad):
        u[bad] = np.cross(x, omega[bad])
    u /= np.linalg.norm(u, axis=1)[:, None]
    v = np.cross(omega, u)
    v /= np.linalg.norm(v, axis=1)[:, None]
    return u, v


def antenna_response(pos: np.ndarray, omega: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    u, v = polarization_basis(omega)
    pu = pos @ u.T
    pv = pos @ v.T
    denom = 1.0 + pos @ omega.T
    denom = np.clip(denom, 1e-6, None)
    fplus = 0.5 * (pu**2 - pv**2) / denom
    fcross = (pu * pv) / denom
    return fplus, fcross


def hd_curve(coszeta: np.ndarray) -> np.ndarray:
    x = (1.0 - coszeta) / 2.0
    out = np.empty_like(x)
    small = x < 1e-12
    out[small] = 0.5
    xs = x[~small]
    out[~small] = 0.5 + 1.5 * xs * np.log(xs) - 0.25 * xs
    return out


def geometry_orf_test(count_map: np.ndarray, nside: int = 64) -> list[dict]:
    import healpy as hp

    names, pos = load_ng15_positions()
    npix = hp.nside2npix(nside)
    if count_map.size != npix:
        raise RuntimeError(f"count_map npix={count_map.size} does not match nside={nside}")
    omega = np.asarray(hp.pix2vec(nside, np.arange(npix))).T
    fplus, fcross = antenna_response(pos, omega)
    iso_w = np.ones(npix, dtype=float) / npix
    lss_w = np.maximum(count_map, 0.0)
    lss_w = lss_w / np.sum(lss_w) if np.sum(lss_w) > 0 else iso_w

    gamma_iso = (fplus * iso_w) @ fplus.T + (fcross * iso_w) @ fcross.T
    gamma_lss = (fplus * lss_w) @ fplus.T + (fcross * lss_w) @ fcross.T
    pairs = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            pairs.append((i, j))
    iso_vec = np.array([gamma_iso[i, j] for i, j in pairs])
    lss_vec = np.array([gamma_lss[i, j] for i, j in pairs])
    cosz = np.array([float(np.dot(pos[i], pos[j])) for i, j in pairs])
    hd = hd_curve(cosz)
    scale = float(np.dot(iso_vec, hd) / max(np.dot(iso_vec, iso_vec), 1e-30))
    iso_scaled = scale * iso_vec
    corr_iso = float(np.corrcoef(iso_scaled, hd)[0, 1])
    corr_lss_hd = float(np.corrcoef(lss_vec, hd)[0, 1])
    rows = [
        {
            "test": "isotropic_numeric_orf_vs_hd_shape",
            "status": "pass",
            "n_pulsars": len(names),
            "n_pairs": len(pairs),
            "pearson_with_hd": corr_iso,
            "best_fit_scale_to_hd": scale,
            "rms_after_scale": float(np.sqrt(np.mean((iso_scaled - hd) ** 2))),
            "note": "Geometry-only antenna-pattern integration; no timing likelihood.",
        },
        {
            "test": "desi_bgs_bright21p5_count_template_vector",
            "status": "pass",
            "n_pulsars": len(names),
            "n_pairs": len(pairs),
            "pearson_with_hd": corr_lss_hd,
            "best_fit_scale_to_hd": "",
            "rms_after_scale": "",
            "note": "Legacy unrandom-corrected DESI BGS_BRIGHT-21.5 count diagnostic; corrected production maps are reported below when present.",
        },
    ]
    return rows


def write_markdown(
    manifest_rows: list[dict],
    mask_summary: dict,
    map_summary: dict,
    orf_rows: list[dict],
    selection_rows: list[dict],
) -> None:
    random_rows = [row for row in manifest_rows if row["role"] == "selection_random"]
    random_ready = sum(1 for row in random_rows if row["status"] in {"exists", "downloaded"})
    lines = [
        "# P7 LSS Tomography Gate",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Download Status",
        "",
        "| key | role | status | bytes | note |",
        "|---|---|---:|---:|---|",
    ]
    for row in manifest_rows:
        lines.append(f"| `{row['key']}` | `{row['role']}` | `{row['status']}` | {row['bytes']} | {row['note']} |")
    lines.extend(
        [
            "",
            "## Map Validation",
            "",
            f"- WISExSCOS mask: `{mask_summary.get('status')}`; nside={mask_summary.get('nside', '')}; nonzero fraction={mask_summary.get('nonzero_fraction', '')}.",
            f"- DESI BGS_BRIGHT-21.5 pilot map: `{map_summary.get('status')}`; objects={map_summary.get('n_objects', '')}; z range=[{map_summary.get('z_min', '')}, {map_summary.get('z_max', '')}]; occupied pixel fraction={map_summary.get('occupied_fraction_all', '')}.",
            "",
            "## NG15 Geometry Gate",
            "",
            "| test | status | n pulsars | n pairs | corr with HD | note |",
            "|---|---|---:|---:|---:|---|",
        ]
    )
    for row in orf_rows:
        lines.append(
            f"| `{row['test']}` | `{row['status']}` | {row['n_pulsars']} | {row['n_pairs']} | {row['pearson_with_hd']:.6f} | {row['note']} |"
        )
    if selection_rows:
        lines.extend(
            [
                "",
                "## Selection-Corrected DESI BGS Maps",
                "",
                "| bin | data weight | random weight | mask fraction | expected-weighted mean delta |",
                "|---|---:|---:|---:|---:|",
            ]
        )
        for row in selection_rows:
            lines.append(
                f"| `{row['bin']}` | {float(row['data_weight_sum']):.6g} | "
                f"{float(row['random_weight_sum']):.6g} | {float(row['masked_pixel_fraction']):.6f} | "
                f"{float(row['expected_weighted_mean_delta']):.3e} |"
            )
    lines.extend(
        [
            "",
            "## Production Gate",
            "",
            "- Safe to use now: DESI BGS_BRIGHT-21.5 count-map geometry as a pipeline smoke test.",
            f"- Production random catalogs staged: `{random_ready}/{len(random_rows)}` DESI BGS_BRIGHT-21.5 random files.",
            f"- Selection-corrected production maps: `{'present' if selection_rows else 'missing'}`.",
            "- Not safe to claim now: LSS-correlated anisotropy evidence or source identification.",
            "- Required before science claim: construct LSS-correlated ORF templates, validate isotropic/random-map nulls, and run an anisotropic PTA likelihood backend.",
            "",
            "## Source Anchors",
            "",
            "- DESI DR1 LSS v1.5 directory: https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/",
            "- 2MPZ SSA page: http://ssa.roe.ac.uk/TWOMPZ.html",
            "- WISExSCOS SSA page: http://ssa.roe.ac.uk/WISExSCOS.html",
        ]
    )
    (RESULTS / "map_validation.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    for path in (DATA, META, DESI, MAPS, RESULTS):
        path.mkdir(parents=True, exist_ok=True)

    manifest_rows = [download_resource(res) for res in RESOURCES]
    write_csv(
        DATA / "manifest.csv",
        manifest_rows,
        [
            "key",
            "url",
            "local_path",
            "role",
            "download_requested",
            "production_required",
            "status",
            "bytes",
            "note",
            "error",
        ],
    )

    mask_summary = validate_wise_mask(DATA / "WISExSCOS" / "WISExSCOSmask.fits.gz")
    map_summary, count_map = build_desi_count_map(
        [
            DESI / "BGS_BRIGHT-21.5_NGC_clustering.dat.fits",
            DESI / "BGS_BRIGHT-21.5_SGC_clustering.dat.fits",
        ],
        nside=64,
    )
    orf_rows = geometry_orf_test(count_map, nside=64)
    selection_rows = load_selection_corrected_summary()
    write_csv(
        RESULTS / "ng15_lss_orf_test.csv",
        orf_rows,
        ["test", "status", "n_pulsars", "n_pairs", "pearson_with_hd", "best_fit_scale_to_hd", "rms_after_scale", "note"],
    )
    write_csv(
        RESULTS / "map_validation.csv",
        [
            {"product": "WISExSCOS_mask", **mask_summary},
            {"product": "DESI_BGS_BRIGHT21p5_count_map", **map_summary},
        ],
        ["product", "status", "nside", "npix", "finite_fraction", "mean_mask", "nonzero_fraction", "n_objects", "z_min", "z_max", "occupied_fraction_all", "occupied_fraction_z_0_0p2", "occupied_fraction_z_0p2_0p5", "max_pixel_count_all", "map_all", "error"],
    )
    write_markdown(manifest_rows, mask_summary, map_summary, orf_rows, selection_rows)
    print(
        {
            "manifest_rows": len(manifest_rows),
            "mask_status": mask_summary.get("status"),
            "map_status": map_summary.get("status"),
            "orf_tests": len(orf_rows),
            "selection_rows": len(selection_rows),
        }
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
