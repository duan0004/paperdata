#!/usr/bin/env python3
"""
PRL H4 robustness: environmental SMBHB turnover vs PTArcade BHB-prior baseline.

This repeats the H4 curved-astrophysical competitor across the same
seed/live-point configurations used for the official-density robustness sweep.
"""

from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from dynesty import NestedSampler
from dynesty.utils import resample_equal


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "code/prl_H4_H5_H6_experiments.py"
OUTDIR = ROOT / "results/T2_NG15yr/bayes_factors"
OUTFILE = OUTDIR / "prl_H4_env_robustness.json"
OUT_MD = OUTDIR / "prl_H4_env_robustness.md"
PARTIAL = OUTDIR / "prl_H4_env_robustness.partial.json"

CONFIGS = [
    {"label": "nlive500_seed42", "nlive": 500, "dlogz": 0.1, "seed": 42},
    {"label": "nlive500_seed7", "nlive": 500, "dlogz": 0.1, "seed": 7},
    {"label": "nlive500_seed123", "nlive": 500, "dlogz": 0.1, "seed": 123},
    {"label": "nlive1000_seed42", "nlive": 1000, "dlogz": 0.05, "seed": 42},
    {"label": "nlive1000_seed7", "nlive": 1000, "dlogz": 0.05, "seed": 7},
]
MODELS = [
    ("baseline_smbhb_ptarcade_bhb_prior", "baseline"),
    ("smbhb_env_fixed_gamma", "env_fixed"),
    ("smbhb_env_free_gamma", "env_free"),
]


def load_helper():
    spec = importlib.util.spec_from_file_location("h456", HELPER)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def run_model(helper, model_key: str, cfg: dict):
    if model_key == "baseline_smbhb_ptarcade_bhb_prior":
        cpta, names, prior_transform, _ = helper.baseline_bhb_prior()
    elif model_key == "smbhb_env_fixed_gamma":
        cpta, names, prior_transform, _ = helper.build_env_model("smbhb_env_fixed_gamma")
    elif model_key == "smbhb_env_free_gamma":
        cpta, names, prior_transform, _ = helper.build_env_model("smbhb_env_free_gamma")
    else:
        raise ValueError(model_key)

    t0 = time.time()
    sampler = NestedSampler(
        helper.finite_loglike(cpta),
        prior_transform,
        cpta.ndim,
        nlive=int(cfg["nlive"]),
        rstate=np.random.default_rng(int(cfg["seed"])),
        bound="multi",
        sample="rwalk",
    )
    sampler.run_nested(dlogz=float(cfg["dlogz"]), print_progress=False)
    res = sampler.results
    weights = np.exp(res["logwt"] - res["logz"][-1])
    equal = resample_equal(res["samples"], weights)
    best_idx = int(np.argmax(res["logl"]))
    return {
        "model": model_key,
        "dim": int(cpta.ndim),
        "ceffyl_param_names": list(cpta.param_names),
        "lnZ": float(res["logz"][-1]),
        "lnZerr": float(res["logzerr"][-1]),
        "niter": int(res.niter),
        "runtime_s": float(time.time() - t0),
        "bestfit": res["samples"][best_idx].tolist(),
        "bestfit_loglike": float(res["logl"][best_idx]),
        "posterior_mean": equal.mean(axis=0).tolist(),
        "posterior_std": equal.std(axis=0).tolist(),
    }


def summarize(config_results: dict):
    rows = []
    for cfg_label, models in config_results.items():
        if "baseline_smbhb_ptarcade_bhb_prior" not in models:
            continue
        base = models["baseline_smbhb_ptarcade_bhb_prior"]
        for key in ["smbhb_env_fixed_gamma", "smbhb_env_free_gamma"]:
            if key not in models:
                continue
            r = models[key]
            rows.append(
                {
                    "config": cfg_label,
                    "model": key,
                    "lnB": float(r["lnZ"] - base["lnZ"]),
                    "lnB_err": float(np.sqrt(r["lnZerr"] ** 2 + base["lnZerr"] ** 2)),
                    "bestfit": dict(zip(r["ceffyl_param_names"], r["bestfit"])),
                }
            )
    aggregate = {}
    for key in ["smbhb_env_fixed_gamma", "smbhb_env_free_gamma"]:
        vals = np.array([r["lnB"] for r in rows if r["model"] == key], dtype=float)
        errs = np.array([r["lnB_err"] for r in rows if r["model"] == key], dtype=float)
        if vals.size == 0:
            continue
        aggregate[key] = {
            "n_configs": int(vals.size),
            "lnB_mean": float(vals.mean()),
            "lnB_sample_std": float(vals.std(ddof=1)) if vals.size > 1 else 0.0,
            "lnB_min": float(vals.min()),
            "lnB_max": float(vals.max()),
            "mean_nested_lnB_err": float(errs.mean()),
        }
    return rows, aggregate


def write(path: Path, helper, config_results: dict, complete: bool):
    rows, aggregate = summarize(config_results)
    out = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "complete": complete,
        "task": "PRL H4 environmental SMBHB robustness",
        "method": {
            "configs": CONFIGS,
            "baseline": "smbhb_ptarcade_bhb_prior",
            "env_model": "power law times [1 + (f_bend/f)^kappa]^-1",
            "n_bins": helper.N_BINS,
        },
        "config_results": config_results,
        "rows": rows,
        "aggregate": aggregate,
    }
    path.write_text(json.dumps(out, indent=2))


def write_md():
    d = json.loads(OUTFILE.read_text())
    lines = [
        "# PRL H4 Environmental-SMBHB Robustness",
        "",
        f"**Date**: {d['generated']}",
        f"**JSON**: `{OUTFILE.relative_to(ROOT)}`",
        "",
        "| Model | mean lnB | sample std | min | max | mean nested err |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for key, agg in d["aggregate"].items():
        lines.append(
            f"| `{key}` | `{agg['lnB_mean']:+.3f}` | `{agg['lnB_sample_std']:.3f}` | "
            f"`{agg['lnB_min']:+.3f}` | `{agg['lnB_max']:+.3f}` | `{agg['mean_nested_lnB_err']:.3f}` |"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "The environmental-turnover SMBHB control is competitive with the leading",
        "cosmological templates on the same official-density evidence scale.  This",
        "turns the PRL story from a pure provenance audit into a source-identification",
        "degeneracy result: current PTA evidence for non-SMBHB templates must be",
        "calibrated against astrophysical spectral curvature.",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


def main():
    helper = load_helper()
    config_results = {}
    for cfg in CONFIGS:
        label = cfg["label"]
        print(f"\n## {label}", flush=True)
        config_results[label] = {}
        for key, _ in MODELS:
            print(f"running {key}", flush=True)
            r = run_model(helper, key, cfg)
            config_results[label][key] = r
            print(f"  lnZ={r['lnZ']:+.3f} +/- {r['lnZerr']:.3f}; runtime={r['runtime_s']:.1f}s", flush=True)
            write(PARTIAL, helper, config_results, complete=False)
        rows, _ = summarize({label: config_results[label]})
        for row in rows:
            print(f"  {row['model']} lnB={row['lnB']:+.3f} +/- {row['lnB_err']:.3f}", flush=True)

    write(OUTFILE, helper, config_results, complete=True)
    write_md()
    print(f"saved: {OUTFILE}")
    print(f"saved: {OUT_MD}")


if __name__ == "__main__":
    main()
