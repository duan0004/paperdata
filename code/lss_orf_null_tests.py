#!/usr/bin/env python3
"""Build LSS-correlated ORF pair templates and random-map null tests.

This script projects DESI BGS_BRIGHT-21.5 selection-corrected overdensity maps
onto NANOGrav 15-year pulsar-pair response vectors.  The outputs are geometry
and null-control products only; they are not timing-residual likelihoods and do
not constitute LSS-correlated anisotropy evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MAP_DIR = ROOT / "data" / "LSS" / "maps"
OUT_DIR = ROOT / "results" / "lss_tomography"
TEMPLATE_DIR = ROOT / "data" / "LSS" / "orf_templates"


@dataclass(frozen=True)
class MapSpec:
    key: str
    delta: Path
    mask: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="prod_nside64_weight_nrand36")
    parser.add_argument("--nside", type=int, default=64)
    parser.add_argument("--nulls", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260423)
    return parser.parse_args()


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def load_ng15_positions() -> tuple[list[str], np.ndarray]:
    import contextlib
    import io

    from enterprise_extensions.load_feathers import load_feathers_from_folder

    feather_dir = ROOT / "data" / "NG15yr" / "tutorials" / "data" / "feathers"
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        psrs = load_feathers_from_folder(str(feather_dir))
    names = []
    pos = []
    for psr in psrs:
        if hasattr(psr, "pos"):
            arr = np.asarray(psr.pos, dtype=np.float64)
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
    bad = np.linalg.norm(u, axis=1) < 1.0e-12
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
    denom = np.clip(denom, 1.0e-6, None)
    fplus = 0.5 * (pu**2 - pv**2) / denom
    fcross = (pu * pv) / denom
    return fplus, fcross


def hd_curve(coszeta: np.ndarray) -> np.ndarray:
    x = (1.0 - coszeta) / 2.0
    out = np.empty_like(x)
    small = x < 1.0e-12
    out[small] = 0.5
    xs = x[~small]
    out[~small] = 0.5 + 1.5 * xs * np.log(xs) - 0.25 * xs
    return out


def pairs_and_hd(names: list[str], pos: np.ndarray) -> tuple[list[tuple[int, int]], np.ndarray, list[dict]]:
    pairs = []
    pair_rows = []
    cosz = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            pairs.append((i, j))
            cz = float(np.dot(pos[i], pos[j]))
            cosz.append(cz)
            pair_rows.append({"i": i, "j": j, "pulsar_i": names[i], "pulsar_j": names[j], "cos_zeta": cz})
    return pairs, hd_curve(np.asarray(cosz, dtype=np.float64)), pair_rows


def map_specs(tag: str) -> list[MapSpec]:
    keys = ["all", "z_0p1_0p2", "z_0p2_0p3", "z_0p3_0p4", "z_0p1_0p4"]
    specs = []
    for key in keys:
        delta = MAP_DIR / f"desi_bgs_bright21p5_{key}_{tag}_delta.npy"
        mask = MAP_DIR / f"desi_bgs_bright21p5_{key}_{tag}_mask.npy"
        specs.append(MapSpec(key, delta, mask))
    return specs


def prepare_field(spec: MapSpec, nside: int) -> tuple[np.ndarray, dict]:
    import healpy as hp

    delta = np.asarray(np.load(spec.delta), dtype=np.float64)
    mask = np.asarray(np.load(spec.mask), dtype=bool)
    npix = hp.nside2npix(nside)
    if delta.size != npix or mask.size != npix:
        raise RuntimeError(f"{spec.key}: npix mismatch for nside={nside}")
    valid = mask & np.isfinite(delta)
    field = np.zeros(npix, dtype=np.float64)
    if np.count_nonzero(valid) == 0:
        raise RuntimeError(f"{spec.key}: empty valid mask")
    field[valid] = delta[valid] - float(np.mean(delta[valid]))
    rms = float(np.sqrt(np.mean(field[valid] ** 2)))
    if rms > 0:
        field[valid] /= rms
    summary = {
        "bin": spec.key,
        "delta_map": rel(spec.delta),
        "mask_map": rel(spec.mask),
        "valid_pixels": int(np.count_nonzero(valid)),
        "mask_fraction": float(np.mean(valid)),
        "area_mean_before_norm": float(np.mean(delta[valid])),
        "area_rms_before_norm": rms,
    }
    return field, summary


def pair_vector(field: np.ndarray, fplus: np.ndarray, fcross: np.ndarray, pairs: list[tuple[int, int]]) -> np.ndarray:
    weight = field / field.size
    gamma = (fplus * weight) @ fplus.T + (fcross * weight) @ fcross.T
    return np.asarray([gamma[i, j] for i, j in pairs], dtype=np.float64)


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size != y.size or x.size < 2:
        return float("nan")
    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def vector_stats(vec: np.ndarray, hd: np.ndarray) -> dict:
    return {
        "pair_mean": float(np.mean(vec)),
        "pair_std": float(np.std(vec)),
        "pair_l2": float(np.linalg.norm(vec)),
        "pearson_with_hd": safe_corr(vec, hd),
        "max_abs_pair": float(np.max(np.abs(vec))),
    }


def null_test(
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
    max_abs = []
    for _ in range(n_nulls):
        shuffled = np.zeros_like(field)
        shuffled_values = values.copy()
        rng.shuffle(shuffled_values)
        shuffled[valid] = shuffled_values
        vec = pair_vector(shuffled, fplus, fcross, pairs)
        stats = vector_stats(vec, hd)
        hd_corrs.append(stats["pearson_with_hd"])
        l2s.append(stats["pair_l2"])
        max_abs.append(stats["max_abs_pair"])
    hd_corrs = np.asarray(hd_corrs, dtype=np.float64)
    l2s = np.asarray(l2s, dtype=np.float64)
    max_abs = np.asarray(max_abs, dtype=np.float64)
    real_corr = float(real_stats["pearson_with_hd"])
    return {
        "n_nulls": n_nulls,
        "null_hd_corr_mean": float(np.nanmean(hd_corrs)),
        "null_hd_corr_std": float(np.nanstd(hd_corrs)),
        "real_hd_corr": real_corr,
        "real_hd_corr_z": float((real_corr - np.nanmean(hd_corrs)) / max(np.nanstd(hd_corrs), 1.0e-30)),
        "real_abs_hd_corr_percentile": float(np.mean(np.abs(hd_corrs) <= abs(real_corr))),
        "null_l2_mean": float(np.mean(l2s)),
        "null_l2_std": float(np.std(l2s)),
        "real_l2_z": float((float(real_stats["pair_l2"]) - np.mean(l2s)) / max(np.std(l2s), 1.0e-30)),
        "null_max_abs_mean": float(np.mean(max_abs)),
        "null_max_abs_std": float(np.std(max_abs)),
    }


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_markdown(tag: str, summary_rows: list[dict], null_rows: list[dict], corr_rows: list[dict], args: argparse.Namespace) -> None:
    lines = [
        "# DESI LSS ORF Template Null Tests",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Scope",
        "",
        "- Input maps: DESI BGS_BRIGHT-21.5 selection-corrected overdensity maps.",
        "- Projection: NANOGrav 15-year pulsar-pair antenna response geometry.",
        "- Claim boundary: geometry/null-control products only; no timing-residual likelihood or source-identification evidence.",
        "",
        "## Configuration",
        "",
        f"- Map tag: `{tag}`.",
        f"- HEALPix nside: `{args.nside}`.",
        f"- Random-map nulls per bin: `{args.nulls}`.",
        f"- RNG seed: `{args.seed}`.",
        "",
        "## Template Summary",
        "",
        "| bin | valid pixels | mask fraction | pair std | corr with HD | ORF vector |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        lines.append(
            f"| `{row['bin']}` | {row['valid_pixels']} | {float(row['mask_fraction']):.6f} | "
            f"{float(row['pair_std']):.6g} | {float(row['pearson_with_hd']):.6f} | `{row['orf_vector']}` |"
        )
    lines.extend(
        [
            "",
            "## Random-Map Nulls",
            "",
            "| bin | real corr(HD) | null mean | null std | z | abs-percentile | real l2 z |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in null_rows:
        lines.append(
            f"| `{row['bin']}` | {float(row['real_hd_corr']):.6f} | {float(row['null_hd_corr_mean']):.6f} | "
            f"{float(row['null_hd_corr_std']):.6f} | {float(row['real_hd_corr_z']):.3f} | "
            f"{float(row['real_abs_hd_corr_percentile']):.3f} | {float(row['real_l2_z']):.3f} |"
        )
    lines.extend(
        [
            "",
            "## Inter-Template Correlations",
            "",
            "| template_a | template_b | pair-vector Pearson r |",
            "|---|---|---:|",
        ]
    )
    for row in corr_rows:
        lines.append(f"| `{row['template_a']}` | `{row['template_b']}` | {float(row['pearson']):.6f} |")
    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "- Passing this test means the DESI maps can be converted into stable pulsar-pair ORF vectors.",
            "- It does not mean that PTA data prefer an LSS-correlated component.",
            "- The next required test is an anisotropic ENTERPRISE likelihood with isotropic and random-map null injections.",
            "",
        ]
    )
    (OUT_DIR / f"lss_orf_null_tests_{tag}.md").write_text("\n".join(lines))


def main() -> int:
    args = parse_args()
    import healpy as hp

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TEMPLATE_DIR.mkdir(parents=True, exist_ok=True)
    names, pos = load_ng15_positions()
    pairs, hd, pair_rows = pairs_and_hd(names, pos)
    npix = hp.nside2npix(args.nside)
    omega = np.asarray(hp.pix2vec(args.nside, np.arange(npix))).T
    fplus, fcross = antenna_response(pos, omega)

    pair_path = TEMPLATE_DIR / f"ng15_pair_index_{args.tag}.csv"
    write_csv(pair_path, pair_rows, ["i", "j", "pulsar_i", "pulsar_j", "cos_zeta"])

    rng = np.random.default_rng(args.seed)
    summary_rows = []
    null_rows = []
    vectors = {}
    for spec in map_specs(args.tag):
        if not spec.delta.exists() or not spec.mask.exists():
            raise FileNotFoundError(f"Missing map inputs for {spec.key}: {spec.delta}, {spec.mask}")
        field, summary = prepare_field(spec, args.nside)
        vec = pair_vector(field, fplus, fcross, pairs)
        stats = vector_stats(vec, hd)
        out_path = TEMPLATE_DIR / f"desi_bgs_bright21p5_{spec.key}_{args.tag}_ng15_orf_pairs.npy"
        np.save(out_path, vec)
        vectors[spec.key] = vec
        row = {**summary, **stats, "orf_vector": rel(out_path), "pair_index": rel(pair_path)}
        summary_rows.append(row)
        null = null_test(field, fplus, fcross, pairs, hd, stats, args.nulls, rng)
        null_rows.append({"bin": spec.key, **null})
        print(f"processed {spec.key}: corr_hd={stats['pearson_with_hd']:.6f}", flush=True)

    corr_rows = []
    keys = list(vectors)
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            corr_rows.append({"template_a": a, "template_b": b, "pearson": safe_corr(vectors[a], vectors[b])})

    summary_fields = [
        "bin",
        "delta_map",
        "mask_map",
        "valid_pixels",
        "mask_fraction",
        "area_mean_before_norm",
        "area_rms_before_norm",
        "pair_mean",
        "pair_std",
        "pair_l2",
        "pearson_with_hd",
        "max_abs_pair",
        "orf_vector",
        "pair_index",
    ]
    write_csv(OUT_DIR / f"lss_orf_template_summary_{args.tag}.csv", summary_rows, summary_fields)
    write_csv(
        OUT_DIR / f"lss_orf_null_summary_{args.tag}.csv",
        null_rows,
        [
            "bin",
            "n_nulls",
            "real_hd_corr",
            "null_hd_corr_mean",
            "null_hd_corr_std",
            "real_hd_corr_z",
            "real_abs_hd_corr_percentile",
            "null_l2_mean",
            "null_l2_std",
            "real_l2_z",
            "null_max_abs_mean",
            "null_max_abs_std",
        ],
    )
    write_csv(OUT_DIR / f"lss_orf_template_correlations_{args.tag}.csv", corr_rows, ["template_a", "template_b", "pearson"])
    write_markdown(args.tag, summary_rows, null_rows, corr_rows, args)
    (OUT_DIR / f"lss_orf_null_run_{args.tag}.json").write_text(
        json.dumps(
            {
                "generated": datetime.now().isoformat(timespec="seconds"),
                "tag": args.tag,
                "nside": args.nside,
                "nulls": args.nulls,
                "seed": args.seed,
                "n_pulsars": len(names),
                "n_pairs": len(pairs),
                "summary": rel(OUT_DIR / f"lss_orf_template_summary_{args.tag}.csv"),
                "null_summary": rel(OUT_DIR / f"lss_orf_null_summary_{args.tag}.csv"),
            },
            indent=2,
        )
        + "\n"
    )
    print({"bins": len(summary_rows), "nulls_per_bin": args.nulls, "n_pairs": len(pairs)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
