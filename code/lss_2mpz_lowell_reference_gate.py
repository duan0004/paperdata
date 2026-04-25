#!/usr/bin/env python3
"""Build fallback low-ell 2MPZ ORF templates for the NG15 null gate.

This script implements the parts of the published NG15+2MPZ LSS-null setup
that are available in the local environment:

* 2MPZ redshift slices 0<z<=0.1 and 0.1<z<=0.2;
* Galactic latitude cut |b|>=20 deg;
* a 2 deg apodization fallback;
* weighted monopole/dipole removal on apodized support;
* low-ell reconstruction with 2<=ell<=12;
* projection onto NANOGrav 15-year pulsar-pair ORF geometry.

Because NaMaster/pymaster is not installed locally, the apodization and
low-ell reconstruction are marked as a fallback implementation.  The outputs
are a production gate for testing the code path; they are not a reproduction of
the published NG15+2MPZ null and must not be interpreted as LSS evidence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
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
class SliceSpec:
    key: str
    z_min: float
    z_max: float

    def contains(self, z: np.ndarray) -> np.ndarray:
        return np.isfinite(z) & (z > self.z_min) & (z <= self.z_max)


SLICES = [
    SliceSpec("z_0_0p1", 0.0, 0.1),
    SliceSpec("z_0p1_0p2", 0.1, 0.2),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nside", type=int, default=64)
    parser.add_argument("--tag", default="prod_nside64_lowell_reference_gate_20260424")
    parser.add_argument("--bcut-deg", type=float, default=20.0)
    parser.add_argument("--apod-fwhm-deg", type=float, default=2.0)
    parser.add_argument("--support-tau", type=float, default=1.0e-3)
    parser.add_argument("--lmax", type=int, default=12)
    parser.add_argument("--nulls", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260424)
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


def load_2mpz_counts(path: Path, nside: int) -> tuple[dict[str, np.ndarray], dict]:
    import healpy as hp

    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)

    maps = {spec.key: np.zeros(hp.nside2npix(nside), dtype=np.float64) for spec in SLICES}
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
            if not (math.isfinite(ra) and math.isfinite(dec) and math.isfinite(z)):
                continue
            if z <= 0.0:
                continue
            rows_valid += 1
            z_min = min(z_min, z)
            z_max = max(z_max, z)
            pix = hp.ang2pix(nside, ra, dec, lonlat=True)
            z_arr = np.asarray([z])
            for spec in SLICES:
                if spec.contains(z_arr)[0]:
                    maps[spec.key][pix] += 1.0
    summary = {
        "catalog": "2MPZ",
        "input": rel(path),
        "rows_read": rows_read,
        "rows_valid": rows_valid,
        "z_min": z_min if math.isfinite(z_min) else "",
        "z_max": z_max if math.isfinite(z_max) else "",
    }
    return maps, summary


def galactic_cut_mask(nside: int, bcut_deg: float) -> np.ndarray:
    import healpy as hp
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    ra = np.degrees(phi)
    dec = 90.0 - np.degrees(theta)
    coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    return np.abs(coord.galactic.b.deg) >= bcut_deg


def apodize_mask(mask: np.ndarray, fwhm_deg: float) -> np.ndarray:
    import healpy as hp

    smooth = hp.smoothing(mask.astype(np.float64), fwhm=np.radians(fwhm_deg), verbose=False)
    smooth = np.asarray(smooth, dtype=np.float64)
    smooth[~np.isfinite(smooth)] = 0.0
    return np.clip(smooth, 0.0, 1.0)


def weighted_overdensity(counts: np.ndarray, weight: np.ndarray, tau: float) -> tuple[np.ndarray, dict]:
    valid = (weight > tau) & np.isfinite(counts)
    delta = np.zeros_like(counts, dtype=np.float64)
    if not np.any(valid):
        raise RuntimeError("empty apodized support")
    total_weight = float(np.sum(weight[valid]))
    mean_count = float(np.sum(weight[valid] * counts[valid]) / total_weight)
    if mean_count <= 0.0:
        raise RuntimeError("non-positive weighted mean count")
    delta[valid] = counts[valid] / mean_count - 1.0
    return delta, {
        "support_pixels": int(np.count_nonzero(valid)),
        "support_fraction": float(np.mean(valid)),
        "weighted_mean_count": mean_count,
        "weighted_delta_mean_before_mono_dipole": float(np.sum(weight[valid] * delta[valid]) / total_weight),
        "delta_rms_before_mono_dipole": float(np.sqrt(np.sum(weight[valid] * delta[valid] ** 2) / total_weight)),
    }


def remove_weighted_monopole_dipole(delta: np.ndarray, weight: np.ndarray, tau: float) -> tuple[np.ndarray, dict]:
    import healpy as hp

    valid = (weight > tau) & np.isfinite(delta)
    pix = np.nonzero(valid)[0]
    if pix.size < 8:
        raise RuntimeError("not enough pixels for monopole/dipole removal")
    nside = hp.npix2nside(delta.size)
    xyz = np.asarray(hp.pix2vec(nside, pix)).T
    design = np.column_stack([np.ones(pix.size), xyz])
    y = delta[pix]
    w = np.sqrt(weight[pix])
    beta, *_ = np.linalg.lstsq(design * w[:, None], y * w, rcond=None)
    clean = np.zeros_like(delta, dtype=np.float64)
    clean[pix] = y - design @ beta
    total_weight = float(np.sum(weight[pix]))
    summary = {
        "monopole_fit": float(beta[0]),
        "dipole_x_fit": float(beta[1]),
        "dipole_y_fit": float(beta[2]),
        "dipole_z_fit": float(beta[3]),
        "weighted_delta_mean_after_mono_dipole": float(np.sum(weight[pix] * clean[pix]) / total_weight),
        "delta_rms_after_mono_dipole": float(np.sqrt(np.sum(weight[pix] * clean[pix] ** 2) / total_weight)),
    }
    return clean, summary


def lowell_fallback(field: np.ndarray, weight: np.ndarray, lmax: int, tau: float) -> tuple[np.ndarray, dict]:
    import healpy as hp

    nside = hp.npix2nside(field.size)
    weighted_field = np.zeros_like(field, dtype=np.float64)
    valid = weight > tau
    weighted_field[valid] = field[valid] * weight[valid]
    alm = hp.map2alm(weighted_field, lmax=lmax, iter=3)
    for ell in (0, 1):
        for emm in range(ell + 1):
            alm[hp.Alm.getidx(lmax, ell, emm)] = 0.0
    lowell = hp.alm2map(alm, nside, lmax=lmax, verbose=False)
    lowell = np.asarray(lowell, dtype=np.float64)
    lowell[~valid] = 0.0
    total_weight = float(np.sum(weight[valid]))
    mean = float(np.sum(weight[valid] * lowell[valid]) / total_weight)
    lowell[valid] -= mean
    rms = float(np.sqrt(np.sum(weight[valid] * lowell[valid] ** 2) / total_weight))
    if rms > 0.0:
        lowell[valid] /= rms
    return lowell, {
        "lowell_lmin": 2,
        "lowell_lmax": lmax,
        "lowell_weighted_mean_before_norm": mean,
        "lowell_weighted_rms_before_norm": rms,
    }


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
    return {
        "n_nulls": n_nulls,
        "real_hd_corr": float(real_stats["pearson_with_hd"]),
        "null_hd_corr_mean": float(np.nanmean(hd_corrs)),
        "null_hd_corr_std": float(np.nanstd(hd_corrs)),
        "real_hd_corr_z": float((float(real_stats["pearson_with_hd"]) - np.nanmean(hd_corrs)) / max(np.nanstd(hd_corrs), 1.0e-30)),
        "real_abs_hd_corr_percentile": float(np.mean(np.abs(hd_corrs) <= abs(float(real_stats["pearson_with_hd"])))),
        "null_l2_mean": float(np.mean(l2s)),
        "null_l2_std": float(np.std(l2s)),
        "real_l2_z": float((float(real_stats["pair_l2"]) - np.mean(l2s)) / max(np.std(l2s), 1.0e-30)),
    }


def build_orf_products(fields: dict[str, np.ndarray], nside: int, tag: str, nulls: int, seed: int) -> tuple[list[dict], list[dict]]:
    import healpy as hp

    ORF_DIR.mkdir(parents=True, exist_ok=True)
    names, pos = lss_orf.load_ng15_positions()
    pairs, hd, pair_rows = lss_orf.pairs_and_hd(names, pos)
    pair_path = ORF_DIR / f"ng15_pair_index_{tag}.csv"
    write_csv(pair_path, pair_rows, ["i", "j", "pulsar_i", "pulsar_j", "cos_zeta"])
    npix = hp.nside2npix(nside)
    omega = np.asarray(hp.pix2vec(nside, np.arange(npix))).T
    fplus, fcross = lss_orf.antenna_response(pos, omega)
    rng = np.random.default_rng(seed)
    summary_rows = []
    null_rows = []
    for key, field in fields.items():
        vec = lss_orf.pair_vector(field, fplus, fcross, pairs)
        stats = lss_orf.vector_stats(vec, hd)
        out_path = ORF_DIR / f"twompz_{key}_{tag}_ng15_orf_pairs.npy"
        np.save(out_path, vec)
        summary_rows.append(
            {
                "bin": key,
                "n_pulsars": len(names),
                "n_pairs": len(pairs),
                **stats,
                "orf_vector": rel(out_path),
                "pair_index": rel(pair_path),
            }
        )
        null_rows.append({"bin": key, **random_map_null(field, fplus, fcross, pairs, hd, stats, nulls, rng)})
    return summary_rows, null_rows


def write_report(args: argparse.Namespace, input_summary: dict, map_rows: list[dict], orf_rows: list[dict], null_rows: list[dict], result: dict) -> None:
    lines = [
        "# 2MPZ Low-Ell Reference-Template Fallback Gate",
        "",
        f"Generated: {result['generated']}",
        "",
        "## Scope",
        "",
        "- Target workflow: NG15+2MPZ published-null reproduction with 2MPZ slices `0<z<=0.1` and `0.1<z<=0.2`.",
        "- Local implementation: Galactic latitude cut, apodized support, weighted monopole/dipole subtraction, low-ell truncation, and NG15 ORF projection.",
        "- Claim boundary: this is a fallback template gate, not a reproduction of the published null and not LSS evidence.",
        "",
        "## Blocking Difference From Reference",
        "",
        "- `pymaster`/NaMaster is not installed locally, so exact C1 apodization and MASTER/MAP low-ell reconstruction are not available.",
        "- This gate uses HEALPix mask smoothing and weighted `map2alm` truncation as a deterministic fallback.",
        "",
        "## Configuration",
        "",
        f"- Tag: `{args.tag}`.",
        f"- HEALPix nside: `{args.nside}`.",
        f"- Galactic cut: `|b| >= {args.bcut_deg} deg`.",
        f"- Apodization fallback FWHM: `{args.apod_fwhm_deg} deg`.",
        f"- Support threshold: `{args.support_tau}`.",
        f"- Low-ell range: `2 <= ell <= {args.lmax}`.",
        f"- Random-map nulls per slice: `{args.nulls}`.",
        "",
        "## 2MPZ Input",
        "",
        "| rows read | rows valid | z min | z max | input |",
        "|---:|---:|---:|---:|---|",
        f"| {input_summary['rows_read']} | {input_summary['rows_valid']} | {input_summary['z_min']} | {input_summary['z_max']} | `{input_summary['input']}` |",
        "",
        "## Map Products",
        "",
        "| slice | counts | support fraction | delta rms | low-ell rms before norm | low-ell map |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in map_rows:
        lines.append(
            f"| `{row['bin']}` | {float(row['counts_sum']):.6g} | {float(row['support_fraction']):.6f} | "
            f"{float(row['delta_rms_after_mono_dipole']):.6g} | {float(row['lowell_weighted_rms_before_norm']):.6g} | `{row['lowell_map']}` |"
        )
    lines.extend(["", "## NG15 Pair ORFs", "", "| slice | corr with HD | pair std | ORF vector |", "|---|---:|---:|---|"])
    for row in orf_rows:
        lines.append(f"| `{row['bin']}` | {float(row['pearson_with_hd']):.6f} | {float(row['pair_std']):.6g} | `{row['orf_vector']}` |")
    lines.extend(["", "## Shuffle Nulls", "", "| slice | real corr(HD) | null mean | null std | z | abs-percentile |", "|---|---:|---:|---:|---:|---:|"])
    for row in null_rows:
        lines.append(
            f"| `{row['bin']}` | {float(row['real_hd_corr']):.6f} | {float(row['null_hd_corr_mean']):.6f} | "
            f"{float(row['null_hd_corr_std']):.6f} | {float(row['real_hd_corr_z']):.3f} | {float(row['real_abs_hd_corr_percentile']):.3f} |"
        )
    lines.extend(
        [
            "",
            "## Decision Boundary",
            "",
            f"- Status: `{result['status']}`.",
            "- The next gate is to run the fixed-parameter NG15 likelihood scan with these ORFs.",
            "- DESI/WISExSCOS remain blocked from science interpretation until the exact NG15+2MPZ published null is reproduced.",
            "",
        ]
    )
    (OUT_DIR / f"lss_2mpz_lowell_reference_gate_{args.tag}.md").write_text("\n".join(lines))


def main() -> int:
    args = parse_args()
    import healpy as hp

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MAP_DIR.mkdir(parents=True, exist_ok=True)
    ORF_DIR.mkdir(parents=True, exist_ok=True)

    twompz_path = DATA / "2MPZ" / "twompz_ra_dec_zphoto.csv.gz"
    counts_by_slice, input_summary = load_2mpz_counts(twompz_path, args.nside)
    binary_mask = galactic_cut_mask(args.nside, args.bcut_deg)
    apod = apodize_mask(binary_mask, args.apod_fwhm_deg)
    np.save(MAP_DIR / f"twompz_galactic_cut_apodized_mask_{args.tag}.npy", apod)

    lowell_fields = {}
    map_rows = []
    for spec in SLICES:
        counts = counts_by_slice[spec.key]
        delta, delta_summary = weighted_overdensity(counts, apod, args.support_tau)
        clean, clean_summary = remove_weighted_monopole_dipole(delta, apod, args.support_tau)
        lowell, lowell_summary = lowell_fallback(clean, apod, args.lmax, args.support_tau)

        stem = MAP_DIR / f"twompz_{spec.key}_{args.tag}"
        counts_path = stem.with_name(stem.name + "_counts.npy")
        delta_path = stem.with_name(stem.name + "_delta.npy")
        clean_path = stem.with_name(stem.name + "_mono_dipole_removed.npy")
        lowell_path = stem.with_name(stem.name + "_lowell.npy")
        np.save(counts_path, counts)
        np.save(delta_path, delta)
        np.save(clean_path, clean)
        np.save(lowell_path, lowell)
        lowell_fields[spec.key] = lowell
        map_rows.append(
            {
                "bin": spec.key,
                "z_min": spec.z_min,
                "z_max": spec.z_max,
                "nside": args.nside,
                "counts_sum": float(np.sum(counts)),
                "binary_mask_fraction": float(np.mean(binary_mask)),
                "apodized_support_fraction": float(np.mean(apod > args.support_tau)),
                "apodized_mask_mean": float(np.mean(apod)),
                **delta_summary,
                **clean_summary,
                **lowell_summary,
                "counts_map": rel(counts_path),
                "delta_map": rel(delta_path),
                "clean_map": rel(clean_path),
                "lowell_map": rel(lowell_path),
                "apodized_mask": rel(MAP_DIR / f"twompz_galactic_cut_apodized_mask_{args.tag}.npy"),
            }
        )
        print(f"built {spec.key}: counts={float(np.sum(counts)):.6g}", flush=True)

    orf_rows, null_rows = build_orf_products(lowell_fields, args.nside, args.tag, args.nulls, args.seed)
    result = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "tag": args.tag,
        "status": "LOWELL_FALLBACK_GATE_PASS",
        "reference_null_reproduction_status": "NOT_REPRODUCED_YET",
        "apodization_method": "healpy_smoothing_fallback_not_exact_C1",
        "lowell_method": "weighted_map2alm_fallback_not_exact_MASTER_MAP",
        "nside": args.nside,
        "lmax": args.lmax,
        "report": rel(OUT_DIR / f"lss_2mpz_lowell_reference_gate_{args.tag}.md"),
    }

    write_csv(OUT_DIR / f"lss_2mpz_lowell_reference_maps_{args.tag}.csv", map_rows, sorted({k for row in map_rows for k in row}))
    write_csv(OUT_DIR / f"lss_2mpz_lowell_reference_orfs_{args.tag}.csv", orf_rows, sorted({k for row in orf_rows for k in row}))
    write_csv(OUT_DIR / f"lss_2mpz_lowell_reference_nulls_{args.tag}.csv", null_rows, sorted({k for row in null_rows for k in row}))
    (OUT_DIR / f"lss_2mpz_lowell_reference_gate_{args.tag}.json").write_text(
        json.dumps(
            {
                "result": result,
                "input": input_summary,
                "maps": map_rows,
                "orfs": orf_rows,
                "nulls": null_rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    write_report(args, input_summary, map_rows, orf_rows, null_rows, result)
    print(json.dumps({"status": result["status"], "report": result["report"]}, indent=2))
    _ = hp  # keep the imported verifier dependency explicit
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
