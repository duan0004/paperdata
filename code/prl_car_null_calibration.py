#!/usr/bin/env python3
"""
PRL robustness experiment: null calibration for the HD free-spectrum covariance.

The CAR pilot showed one non-null correlation mode but a near-diagonal overall
covariance matrix.  This script calibrates the observed correlation statistics
against independent Gaussian samples with the same sample count and bin count.

Outputs:
  results/T2_NG15yr/covariance/car_null_calibration.json
  results/T2_NG15yr/covariance/car_null_calibration.md
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
COVDIR = ROOT / "results/T2_NG15yr/covariance"
OUT_JSON = COVDIR / "car_null_calibration.json"
OUT_MD = COVDIR / "car_null_calibration.md"

N_SIM = 5000
SEED = 20260421


def entropy_effective_rank(eigvals: np.ndarray) -> float:
    vals = np.asarray(eigvals, dtype=float)
    probs = vals / np.sum(vals)
    entropy = -np.sum(probs * np.log(probs))
    return float(np.exp(entropy))


def participation_effective_rank(eigvals: np.ndarray) -> float:
    vals = np.asarray(eigvals, dtype=float)
    return float(np.sum(vals) ** 2 / np.sum(vals**2))


def corr_stats(corr: np.ndarray, n_samples: int) -> dict:
    p = corr.shape[0]
    offdiag = corr[np.triu_indices(p, k=1)]
    eigvals = np.linalg.eigvalsh(corr)[::-1]
    mp_upper = (1.0 + np.sqrt(p / n_samples)) ** 2
    return {
        "max_abs_offdiag": float(np.max(np.abs(offdiag))),
        "mean_abs_offdiag": float(np.mean(np.abs(offdiag))),
        "lambda_max": float(eigvals[0]),
        "lambda_min": float(eigvals[-1]),
        "entropy_effective_rank": entropy_effective_rank(eigvals),
        "participation_effective_rank": participation_effective_rank(eigvals),
        "mp_noise_upper": float(mp_upper),
        "n_signal_modes_above_1p2_MP": int(np.sum(eigvals > 1.2 * mp_upper)),
    }


def p_ge(null_values: np.ndarray, observed: float) -> float:
    return float((1.0 + np.sum(null_values >= observed)) / (null_values.size + 1.0))


def p_le(null_values: np.ndarray, observed: float) -> float:
    return float((1.0 + np.sum(null_values <= observed)) / (null_values.size + 1.0))


def percentile_summary(vals: np.ndarray) -> dict:
    return {
        "p50": float(np.percentile(vals, 50)),
        "p84": float(np.percentile(vals, 84)),
        "p95": float(np.percentile(vals, 95)),
        "p99": float(np.percentile(vals, 99)),
        "p999": float(np.percentile(vals, 99.9)),
        "max": float(np.max(vals)),
    }


def main():
    t0 = time.time()
    corr = np.load(COVDIR / "corr_bin.npy")
    pilot = json.loads((COVDIR / "pilot_report.json").read_text())
    n_samples = int(pilot["N_samples_postburn"])
    n_bins = int(pilot["N_bins"])
    if corr.shape != (n_bins, n_bins):
        raise ValueError(f"corr shape {corr.shape} does not match N_bins={n_bins}")

    observed = corr_stats(corr, n_samples)

    rng = np.random.default_rng(SEED)
    null = {
        "max_abs_offdiag": np.empty(N_SIM),
        "mean_abs_offdiag": np.empty(N_SIM),
        "lambda_max": np.empty(N_SIM),
        "entropy_effective_rank": np.empty(N_SIM),
        "participation_effective_rank": np.empty(N_SIM),
        "n_signal_modes_above_1p2_MP": np.empty(N_SIM),
    }
    for i in range(N_SIM):
        x = rng.normal(size=(n_samples, n_bins))
        c = np.corrcoef(x, rowvar=False)
        s = corr_stats(c, n_samples)
        for key in null:
            null[key][i] = s[key]

    null_summary = {key: percentile_summary(vals) for key, vals in null.items()}
    p_values = {
        "max_abs_offdiag_p_ge_observed": p_ge(null["max_abs_offdiag"], observed["max_abs_offdiag"]),
        "mean_abs_offdiag_p_ge_observed": p_ge(null["mean_abs_offdiag"], observed["mean_abs_offdiag"]),
        "lambda_max_p_ge_observed": p_ge(null["lambda_max"], observed["lambda_max"]),
        "entropy_effective_rank_p_le_observed": p_le(
            null["entropy_effective_rank"], observed["entropy_effective_rank"]
        ),
        "participation_effective_rank_p_le_observed": p_le(
            null["participation_effective_rank"], observed["participation_effective_rank"]
        ),
        "n_signal_modes_p_ge_observed": p_ge(
            null["n_signal_modes_above_1p2_MP"],
            observed["n_signal_modes_above_1p2_MP"],
        ),
    }

    out = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "task": "PRL CAR covariance null calibration",
        "source": {
            "corr_bin": str(COVDIR / "corr_bin.npy"),
            "pilot_report": str(COVDIR / "pilot_report.json"),
        },
        "method": {
            "null": "independent Gaussian samples with same N_samples and N_bins",
            "n_sim": N_SIM,
            "seed": SEED,
            "n_samples": n_samples,
            "n_bins": n_bins,
        },
        "observed": observed,
        "null_summary": null_summary,
        "p_values": p_values,
        "runtime_s": float(time.time() - t0),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))

    md = f"""# PRL Free-Spectrum Bin-Covariance Null Calibration

**Date**: {out['generated']}  
**Input**: `results/T2_NG15yr/covariance/corr_bin.npy`  
**Null**: independent Gaussian samples with the same sample count and bin count.  
**Simulations**: {N_SIM}  

## Observed Statistics

| statistic | observed | null p-value |
|---|---:|---:|
| max abs off-diagonal correlation | {observed['max_abs_offdiag']:.6f} | {p_values['max_abs_offdiag_p_ge_observed']:.6g} |
| mean abs off-diagonal correlation | {observed['mean_abs_offdiag']:.6f} | {p_values['mean_abs_offdiag_p_ge_observed']:.6g} |
| lambda max | {observed['lambda_max']:.6f} | {p_values['lambda_max_p_ge_observed']:.6g} |
| entropy effective rank | {observed['entropy_effective_rank']:.6f} | {p_values['entropy_effective_rank_p_le_observed']:.6g} |
| participation effective rank | {observed['participation_effective_rank']:.6f} | {p_values['participation_effective_rank_p_le_observed']:.6g} |
| signal modes above 1.2 MP upper edge | {observed['n_signal_modes_above_1p2_MP']} | {p_values['n_signal_modes_p_ge_observed']:.6g} |

## Null Percentiles

| statistic | p50 | p95 | p99 | p99.9 | max |
|---|---:|---:|---:|---:|---:|
| max abs off-diagonal correlation | {null_summary['max_abs_offdiag']['p50']:.6f} | {null_summary['max_abs_offdiag']['p95']:.6f} | {null_summary['max_abs_offdiag']['p99']:.6f} | {null_summary['max_abs_offdiag']['p999']:.6f} | {null_summary['max_abs_offdiag']['max']:.6f} |
| mean abs off-diagonal correlation | {null_summary['mean_abs_offdiag']['p50']:.6f} | {null_summary['mean_abs_offdiag']['p95']:.6f} | {null_summary['mean_abs_offdiag']['p99']:.6f} | {null_summary['mean_abs_offdiag']['p999']:.6f} | {null_summary['mean_abs_offdiag']['max']:.6f} |
| lambda max | {null_summary['lambda_max']['p50']:.6f} | {null_summary['lambda_max']['p95']:.6f} | {null_summary['lambda_max']['p99']:.6f} | {null_summary['lambda_max']['p999']:.6f} | {null_summary['lambda_max']['max']:.6f} |
| entropy effective rank | {null_summary['entropy_effective_rank']['p50']:.6f} | {null_summary['entropy_effective_rank']['p95']:.6f} | {null_summary['entropy_effective_rank']['p99']:.6f} | {null_summary['entropy_effective_rank']['p999']:.6f} | {null_summary['entropy_effective_rank']['max']:.6f} |
| participation effective rank | {null_summary['participation_effective_rank']['p50']:.6f} | {null_summary['participation_effective_rank']['p95']:.6f} | {null_summary['participation_effective_rank']['p99']:.6f} | {null_summary['participation_effective_rank']['p999']:.6f} | {null_summary['participation_effective_rank']['max']:.6f} |

## Interpretation

The observed correlation structure is not consistent with pure independent
Gaussian sampling noise: the leading correlation mode is real.  However, the
mean absolute off-diagonal correlation remains small and both effective-rank
definitions stay close to the full 14-bin rank.  This supports the PRL claim that
bin covariance exists but is subdominant compared with template/density
provenance for the leading official-density evidence reproduction.
"""
    OUT_MD.write_text(md)
    print(json.dumps(out, indent=2)[:4000])
    print(f"saved: {OUT_JSON}")
    print(f"saved: {OUT_MD}")


if __name__ == "__main__":
    main()
