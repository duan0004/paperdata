#!/usr/bin/env python3
"""Build NaMaster/MAP 2MPZ low-ell maps for the NG15+2MPZ null gate.

This is a map-construction production gate, not a timing-likelihood or Bayes
factor calculation.  It replaces the earlier HEALPix smoothing fallback with
NaMaster C1 apodization, a MASTER target spectrum, a low-ell MAP/Wiener-style
linear reconstruction, and per-ell target-spectrum rescaling.

The generated maps still require projection into PTA ORF pair vectors and then
an NG15 likelihood/hypermodel run before they can be compared to the published
NG15+2MPZ null.
"""

from __future__ import annotations

import argparse
import csv
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
OUT_DIR = ROOT / "results" / "lss_tomography"

if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import lss_2mpz_lowell_reference_gate as fallback  # noqa: E402


@dataclass(frozen=True)
class Basis:
    labels: list[tuple[int, int, str]]
    matrix: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nside", type=int, default=64)
    parser.add_argument("--tag", default="prod_nside64_namaster_reference_gate_20260424")
    parser.add_argument("--bcut-deg", type=float, default=20.0)
    parser.add_argument("--apod-deg", type=float, default=2.0)
    parser.add_argument("--support-tau", type=float, default=1.0e-3)
    parser.add_argument("--lmax", type=int, default=12)
    parser.add_argument("--prior-floor", type=float, default=1.0e-12)
    parser.add_argument("--ridge", type=float, default=1.0e-10)
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


def namaster_apodize(mask: np.ndarray, apod_deg: float) -> np.ndarray:
    import pymaster as nmt

    apod = nmt.mask_apodization(mask.astype(np.float64), apod_deg, apotype="C1")
    apod = np.asarray(apod, dtype=np.float64)
    apod[~np.isfinite(apod)] = 0.0
    return np.clip(apod, 0.0, 1.0)


def master_target_cl(clean: np.ndarray, apod: np.ndarray, lmax: int, prior_floor: float) -> tuple[np.ndarray, dict]:
    import healpy as hp
    import pymaster as nmt

    nside = hp.npix2nside(clean.size)
    field = nmt.NmtField(apod, [clean])
    bins = nmt.NmtBin.from_nside_linear(nside, 1)
    workspace = nmt.NmtWorkspace()
    workspace.compute_coupling_matrix(field, field, bins)
    coupled = nmt.compute_coupled_cell(field, field)
    decoupled = workspace.decouple_cell(coupled)[0]
    ell_eff = np.asarray(bins.get_effective_ells(), dtype=np.float64)
    ells = np.arange(lmax + 1, dtype=np.float64)
    target = np.interp(ells, ell_eff, decoupled, left=decoupled[0], right=decoupled[-1])
    positive = target[np.isfinite(target) & (target > 0.0)]
    floor = max(prior_floor, float(np.nanmedian(positive)) * 1.0e-6 if positive.size else prior_floor)
    target = np.where(np.isfinite(target) & (target > floor), target, floor)
    target[0:2] = floor
    return target, {
        "target_cl_floor": floor,
        "target_cl_min_l2_lmax": float(np.min(target[2:])),
        "target_cl_max_l2_lmax": float(np.max(target[2:])),
        "namaster_decoupled_bands": int(decoupled.size),
    }


def real_sph_basis(nside: int, lmax: int, pixels: np.ndarray | None = None) -> Basis:
    import healpy as hp
    from scipy.special import sph_harm

    if pixels is None:
        pixels = np.arange(hp.nside2npix(nside))
    theta, phi = hp.pix2ang(nside, pixels)
    cols = []
    labels: list[tuple[int, int, str]] = []
    sqrt2 = math.sqrt(2.0)
    for ell in range(2, lmax + 1):
        y0 = sph_harm(0, ell, phi, theta).real
        cols.append(y0)
        labels.append((ell, 0, "real"))
        for emm in range(1, ell + 1):
            y = sph_harm(emm, ell, phi, theta)
            cols.append(sqrt2 * y.real)
            labels.append((ell, emm, "real"))
            cols.append(sqrt2 * y.imag)
            labels.append((ell, emm, "imag"))
    return Basis(labels=labels, matrix=np.column_stack(cols).astype(np.float64))


def solve_map_coefficients(
    clean: np.ndarray,
    apod: np.ndarray,
    lmax: int,
    target_cl: np.ndarray,
    support_tau: float,
    ridge: float,
) -> tuple[np.ndarray, Basis, dict]:
    import healpy as hp

    nside = hp.npix2nside(clean.size)
    valid = (apod > support_tau) & np.isfinite(clean)
    pix = np.nonzero(valid)[0]
    basis = real_sph_basis(nside, lmax, pix)
    y = clean[pix]
    w = apod[pix]
    prior_var = np.asarray([target_cl[ell] for ell, _, _ in basis.labels], dtype=np.float64)
    prior_prec = 1.0 / np.maximum(prior_var, ridge)
    yw = basis.matrix * np.sqrt(w)[:, None]
    lhs = yw.T @ yw + np.diag(prior_prec + ridge)
    rhs = basis.matrix.T @ (w * y)
    cond = float(np.linalg.cond(lhs))
    coeff = np.linalg.solve(lhs, rhs)
    return coeff, basis, {
        "map_valid_pixels": int(pix.size),
        "map_basis_size": int(len(basis.labels)),
        "map_linear_system_condition": cond,
    }


def coefficients_to_map(nside: int, coeff: np.ndarray, labels: list[tuple[int, int, str]]) -> np.ndarray:
    basis_all = real_sph_basis(nside, max(ell for ell, _, _ in labels))
    if basis_all.labels != labels:
        raise RuntimeError("basis label mismatch")
    return basis_all.matrix @ coeff


def per_ell_rescale(
    nside: int,
    coeff: np.ndarray,
    labels: list[tuple[int, int, str]],
    target_cl: np.ndarray,
    lmax: int,
) -> tuple[np.ndarray, dict]:
    import healpy as hp

    lowell = coefficients_to_map(nside, coeff, labels)
    cl = hp.anafast(lowell, lmax=lmax)
    scaled = coeff.copy()
    factors = {}
    for ell in range(2, lmax + 1):
        if cl[ell] > 0.0 and target_cl[ell] > 0.0:
            factor = math.sqrt(float(target_cl[ell] / cl[ell]))
        else:
            factor = 1.0
        factors[str(ell)] = factor
        for idx, (label_ell, _, _) in enumerate(labels):
            if label_ell == ell:
                scaled[idx] *= factor
    rescaled_map = coefficients_to_map(nside, scaled, labels)
    cl_after = hp.anafast(rescaled_map, lmax=lmax)
    return scaled, {
        "rescale_factors": factors,
        "cl_before_l2_lmax": [float(x) for x in cl[2 : lmax + 1]],
        "cl_after_l2_lmax": [float(x) for x in cl_after[2 : lmax + 1]],
    }


def normalize_on_support(field: np.ndarray, weight: np.ndarray, tau: float) -> tuple[np.ndarray, dict]:
    valid = (weight > tau) & np.isfinite(field)
    out = np.zeros_like(field, dtype=np.float64)
    total_weight = float(np.sum(weight[valid]))
    mean = float(np.sum(weight[valid] * field[valid]) / total_weight)
    centered = field[valid] - mean
    rms = float(np.sqrt(np.sum(weight[valid] * centered**2) / total_weight))
    if rms > 0.0:
        out[valid] = centered / rms
    else:
        out[valid] = centered
    return out, {
        "normalized_weighted_mean_before": mean,
        "normalized_weighted_rms_before": rms,
    }


def write_report(args: argparse.Namespace, input_summary: dict, rows: list[dict], result: dict) -> None:
    lines = [
        "# 2MPZ NaMaster/MAP Map Gate",
        "",
        f"Generated: {result['generated']}",
        "",
        "## Scope",
        "",
        "- Target workflow: NG15+2MPZ published-null map construction for the slices `0<z<=0.1` and `0.1<z<=0.2`.",
        "- This run uses NaMaster C1 apodization, MASTER target spectra, low-ell MAP/Wiener reconstruction, and per-ell target-spectrum rescaling.",
        "- Claim boundary: this is map construction only. It is not a PTA likelihood, not a Bayes factor, and not an LSS detection.",
        "",
        "## Configuration",
        "",
        f"- Tag: `{args.tag}`.",
        f"- HEALPix nside: `{args.nside}`.",
        f"- Galactic cut: `|b| >= {args.bcut_deg} deg`.",
        f"- C1 apodization: `{args.apod_deg} deg`.",
        f"- Support threshold: `{args.support_tau}`.",
        f"- Low-ell range: `2 <= ell <= {args.lmax}`.",
        "",
        "## 2MPZ Input",
        "",
        "| rows read | rows valid | z min | z max | input |",
        "|---:|---:|---:|---:|---|",
        f"| {input_summary['rows_read']} | {input_summary['rows_valid']} | {input_summary['z_min']} | {input_summary['z_max']} | `{input_summary['input']}` |",
        "",
        "## Map Products",
        "",
        "| slice | counts | support fraction | basis size | cond(A) | target Cl min | target Cl max | output map |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['bin']}` | {float(row['counts_sum']):.6g} | {float(row['support_fraction']):.6f} | "
            f"{int(row['map_basis_size'])} | {float(row['map_linear_system_condition']):.6g} | "
            f"{float(row['target_cl_min_l2_lmax']):.6g} | {float(row['target_cl_max_l2_lmax']):.6g} | `{row['namaster_map']}` |"
        )
    lines.extend(
        [
            "",
            "## Decision Boundary",
            "",
            f"- Status: `{result['status']}`.",
            "- The next gate is to project these maps into NG15 pair ORF vectors and rerun the NG15+2MPZ null-reproduction likelihood gate.",
            "- DESI/WISExSCOS remain blocked from science interpretation until the published 2MPZ null is reproduced.",
            "",
        ]
    )
    (OUT_DIR / f"lss_2mpz_namaster_map_gate_{args.tag}.md").write_text("\n".join(lines))


def main() -> int:
    args = parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MAP_DIR.mkdir(parents=True, exist_ok=True)

    twompz_path = DATA / "2MPZ" / "twompz_ra_dec_zphoto.csv.gz"
    counts_by_slice, input_summary = fallback.load_2mpz_counts(twompz_path, args.nside)
    binary_mask = fallback.galactic_cut_mask(args.nside, args.bcut_deg)
    apod = namaster_apodize(binary_mask.astype(np.float64), args.apod_deg)
    mask_path = MAP_DIR / f"twompz_galactic_cut_c1_apodized_mask_{args.tag}.npy"
    np.save(mask_path, apod)

    rows: list[dict] = []
    for spec in fallback.SLICES:
        counts = counts_by_slice[spec.key]
        delta, delta_summary = fallback.weighted_overdensity(counts, apod, args.support_tau)
        clean, clean_summary = fallback.remove_weighted_monopole_dipole(delta, apod, args.support_tau)
        target_cl, target_summary = master_target_cl(clean, apod, args.lmax, args.prior_floor)
        coeff, basis, solve_summary = solve_map_coefficients(
            clean, apod, args.lmax, target_cl, args.support_tau, args.ridge
        )
        scaled_coeff, rescale_summary = per_ell_rescale(args.nside, coeff, basis.labels, target_cl, args.lmax)
        reconstructed = coefficients_to_map(args.nside, scaled_coeff, basis.labels)
        normalized, norm_summary = normalize_on_support(reconstructed, apod, args.support_tau)

        stem = MAP_DIR / f"twompz_{spec.key}_{args.tag}"
        counts_path = stem.with_name(stem.name + "_counts.npy")
        delta_path = stem.with_name(stem.name + "_delta.npy")
        clean_path = stem.with_name(stem.name + "_mono_dipole_removed.npy")
        target_path = stem.with_name(stem.name + "_target_cl.npy")
        coeff_path = stem.with_name(stem.name + "_map_coefficients.npy")
        lowell_path = stem.with_name(stem.name + "_namaster_map.npy")
        np.save(counts_path, counts)
        np.save(delta_path, delta)
        np.save(clean_path, clean)
        np.save(target_path, target_cl)
        np.save(coeff_path, scaled_coeff)
        np.save(lowell_path, normalized)

        rows.append(
            {
                "bin": spec.key,
                "z_min": spec.z_min,
                "z_max": spec.z_max,
                "nside": args.nside,
                "counts_sum": float(np.sum(counts)),
                "binary_mask_fraction": float(np.mean(binary_mask)),
                "apodized_support_fraction": float(np.mean(apod > args.support_tau)),
                "apodized_mask_mean": float(np.mean(apod)),
                "apodized_mask": rel(mask_path),
                **delta_summary,
                **clean_summary,
                **target_summary,
                **solve_summary,
                **norm_summary,
                "rescale_factors_json": json.dumps(rescale_summary["rescale_factors"], sort_keys=True),
                "counts_map": rel(counts_path),
                "delta_map": rel(delta_path),
                "clean_map": rel(clean_path),
                "target_cl": rel(target_path),
                "map_coefficients": rel(coeff_path),
                "namaster_map": rel(lowell_path),
            }
        )
        print(f"built {spec.key}: counts={float(np.sum(counts)):.6g}", flush=True)

    result = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "tag": args.tag,
        "status": "NAMASTER_MAP_GATE_PASS",
        "reference_null_reproduction_status": "NOT_REPRODUCED_YET",
        "nside": args.nside,
        "lmax": args.lmax,
        "report": rel(OUT_DIR / f"lss_2mpz_namaster_map_gate_{args.tag}.md"),
    }
    fields = sorted({key for row in rows for key in row})
    write_csv(OUT_DIR / f"lss_2mpz_namaster_map_gate_{args.tag}.csv", rows, fields)
    (OUT_DIR / f"lss_2mpz_namaster_map_gate_{args.tag}.json").write_text(
        json.dumps({"result": result, "input": input_summary, "maps": rows}, indent=2, sort_keys=True) + "\n"
    )
    write_report(args, input_summary, rows, result)
    print(json.dumps({"status": result["status"], "report": result["report"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

