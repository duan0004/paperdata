#!/usr/bin/env python3
"""Check whether reported official-density posterior samples hit ceffyl grid plateaus.

ceffyl 1.41.2 clips spectra below the lower rho-grid boundary to the bottom
bin and returns -inf if any modeled rho is above the grid.  This diagnostic
counts those events for the posterior samples used in the reported
official-density rows.
"""

from __future__ import annotations

import importlib.metadata as md
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import ptarcade
from ceffyl import Ceffyl
from enterprise.signals import parameter
from enterprise.signals.parameter import function


ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "results/T2_NG15yr/bayes_factors/prl_ceffyl_plateau_diagnostic.json"
OUT_MD = ROOT / "results/T2_NG15yr/bayes_factors/prl_ceffyl_plateau_diagnostic.md"

BRIDGE_PATH = ROOT / "code/prl_reference_bridge_pipeline.py"
spec = importlib.util.spec_from_file_location("bridge_pipeline_for_plateau", BRIDGE_PATH)
bridge = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = bridge
assert spec.loader is not None
spec.loader.exec_module(bridge)

TOP_ROWS = [
    "SIGW-Gaussian",
    "SIGW-Delta",
    "Cosmic-Superstrings",
    "SMBHB-Eccentric",
    "SMBHB-BrokenPL",
    "SMBHB-Env",
]


def make_cpta_signal(model, n_freq: int = 14):
    cpta = Ceffyl.ceffyl(str(bridge.DATADIR_OFFICIAL))

    @function
    def psd(f, Tspan, **kwargs):
        params = dict(model.param_defaults)
        for key in model.param_names:
            if key in kwargs:
                params[key] = kwargs[key]
        val = model.timing_residual_psd(np.asarray(f, dtype=float), **params) / Tspan
        return np.asarray(val, dtype=float)

    params = []
    for key in model.param_names:
        if isinstance(model.param_defaults.get(key), str):
            continue
        lo, hi = model.param_priors[key]
        params.append(parameter.Uniform(float(lo), float(hi))(key))

    sig = Ceffyl.signal(N_freqs=n_freq, psd=psd, params=params, name="")
    cpta.add_signals([sig])
    return cpta


def grid_counts(cpta, xs: np.ndarray) -> tuple[int, int, int, float, float]:
    red_rho = np.zeros((cpta.N_psrs, cpta.N_freqs))
    cp_rho = np.zeros((cpta.N_psrs, cpta.N_freqs))
    for s in cpta.cp_signals:
        mapped_xs = {s_i.name: xs[p] for p, s_i in zip(s.pmap, s.params)}
        cp_rho[s.ixgrid] += s.get_rho(
            cpta.reshaped_freqs[s.freq_idxs],
            Tspan=cpta.Tspan,
            mapped_xs=mapped_xs,
        )
    rho = red_rho + cp_rho
    logrho = 0.5 * np.log10(rho)
    idx = np.searchsorted(cpta.binedges, logrho) - 1
    return (
        int((idx < 0).sum()),
        int((idx >= cpta.rho_grid.shape[0]).sum()),
        int(idx.size),
        float(np.nanmin(logrho)),
        float(np.nanmax(logrho)),
    )


def diagnose_model(model_key: str, model) -> dict:
    cpta = make_cpta_signal(model, n_freq=14)
    order = list(cpta.param_names)
    sample_path = ROOT / "results/prl_reference_bridge/P1_ng15_official" / model_key / "posterior_equal_weight.npz"
    data = np.load(sample_path, allow_pickle=True)
    samples = data["samples"]
    names = [str(x) for x in data["param_names"]]
    indices = np.linspace(0, samples.shape[0] - 1, min(512, samples.shape[0]), dtype=int)
    lower = upper = total = 0
    min_logrho = []
    max_logrho = []
    for row in samples[indices]:
        params = dict(model.param_defaults)
        for name, val in zip(names, row):
            params[name] = float(val)
        xs = np.asarray([params[k] for k in order], dtype=float)
        lo, hi, n_cell, mn, mx = grid_counts(cpta, xs)
        lower += lo
        upper += hi
        total += n_cell
        min_logrho.append(mn)
        max_logrho.append(mx)
    return {
        "model": model_key,
        "n_posterior_samples_checked": int(len(indices)),
        "n_grid_cells_checked": int(total),
        "lower_boundary_fraction": float(lower / total),
        "upper_boundary_fraction": float(upper / total),
        "min_log10_rho": float(min(min_logrho)),
        "max_log10_rho": float(max(max_logrho)),
        "grid_min_log10_rho": float(cpta.binedges[0]),
        "grid_max_log10_rho": float(cpta.binedges[-1]),
    }


def main() -> None:
    models = bridge.bridge_models()
    rows = [diagnose_model(key, models[key]) for key in TOP_ROWS]
    out = {
        "generated": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
        "ptarcade_version": ptarcade.__version__,
        "ceffyl_version": md.version("ceffyl"),
        "official_model_archive": "PTArcade_models_1.0.0/models_1.0.0",
        "official_density_archive": "PTArcade_ceffyl_0.2.0/ng15_30f_fs{hd}_ceffyl",
        "diagnostic": "counts lower-grid clipping and upper-grid -inf events for equal-weight official-density posterior samples",
        "rows": rows,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=2) + "\n")

    md_lines = [
        "# PRL ceffyl plateau diagnostic",
        "",
        f"Generated: {out['generated']}",
        "",
        f"- PTArcade version: `{out['ptarcade_version']}`",
        f"- ceffyl version: `{out['ceffyl_version']}`",
        f"- Official density archive: `{out['official_density_archive']}`",
        "",
        "| Model | samples | lower-boundary fraction | upper-boundary fraction | log10 rho range | grid range |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        md_lines.append(
            f"| `{row['model']}` | {row['n_posterior_samples_checked']} | "
            f"{row['lower_boundary_fraction']:.6f} | {row['upper_boundary_fraction']:.6f} | "
            f"[{row['min_log10_rho']:.3f}, {row['max_log10_rho']:.3f}] | "
            f"[{row['grid_min_log10_rho']:.3f}, {row['grid_max_log10_rho']:.3f}] |"
        )
    OUT_MD.write_text("\n".join(md_lines) + "\n")
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
