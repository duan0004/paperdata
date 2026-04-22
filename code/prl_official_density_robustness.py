#!/usr/bin/env python3
"""
PRL robustness experiment: official PTArcade template + official ceffyl density.

This script repeats the decisive PRL evidence check across several dynesty
settings.  It is intentionally narrow: the goal is to test whether the
official-density reproduction of the NANOGrav 15-year leading SIGW/string
Bayes factors is stable against sampler seed and live-point count.

Outputs:
  results/T2_NG15yr/bayes_factors/prl_official_density_robustness.json

The calculation uses the official PTArcade ceffyl density product from
Zenodo DOI 10.5281/zenodo.10495907 and the official NANOGrav/PTArcade model
files from DOI 10.5281/zenodo.8084351.
"""

from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from ceffyl import Ceffyl, models
from dynesty import NestedSampler
from dynesty.utils import resample_equal
from enterprise.signals import parameter
from scipy.special import ndtri

import ptarcade
import ptarcade.models_utils as aux
from ptarcade.models_utils import ParamDict
from ptarcade.signal_builder import bhb_priors, powerlaw2_ceffyl


ROOT = Path(__file__).resolve().parents[1]
DATADIR = ROOT / "data/NG15yr/PTArcade_ceffyl_0.2.0/ng15_30f_fs{hd}_ceffyl"
MODEL_ROOT = ROOT / "data/NG15yr/PTArcade_models_1.0.0/models_1.0.0"
OUTDIR = ROOT / "results/T2_NG15yr/bayes_factors"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTFILE = OUTDIR / "prl_official_density_robustness.json"
PARTIAL = OUTDIR / "prl_official_density_robustness.partial.json"

N_BINS = 14
PUBLISHED_LNB = {
    "sigw_gauss_ptarcade": float(np.log(57.0)),
    "sigw_delta_ptarcade": float(np.log(44.0)),
    "super_ptarcade": float(np.log(46.0)),
}
PUBLISHED_B = {
    "sigw_gauss_ptarcade": 57.0,
    "sigw_delta_ptarcade": 44.0,
    "super_ptarcade": 46.0,
}

MODEL_FILES = {
    "sigw_gauss_ptarcade": "sigw_gauss.py",
    "sigw_delta_ptarcade": "sigw_delta.py",
    "super_ptarcade": "super.py",
}

CONFIGS = [
    {"label": "nlive500_seed42", "nlive": 500, "dlogz": 0.1, "seed": 42},
    {"label": "nlive500_seed7", "nlive": 500, "dlogz": 0.1, "seed": 7},
    {"label": "nlive500_seed123", "nlive": 500, "dlogz": 0.1, "seed": 123},
    {"label": "nlive1000_seed42", "nlive": 1000, "dlogz": 0.05, "seed": 42},
    {"label": "nlive1000_seed7", "nlive": 1000, "dlogz": 0.05, "seed": 7},
]

KINDS = [
    "smbhb_ptarcade_bhb_prior",
    "sigw_gauss_ptarcade",
    "sigw_delta_ptarcade",
    "super_ptarcade",
]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def build_ceffyl_model(kind: str):
    cpta = Ceffyl.ceffyl(str(DATADIR))
    prior_transform = None
    model_path = None

    if kind == "smbhb_fixed_gamma":
        sig = Ceffyl.signal(
            N_freqs=N_BINS,
            psd=models.powerlaw,
            params=[parameter.Uniform(-18, -11)("log10_A")],
            const_params={"gamma": 13.0 / 3.0},
            name="",
        )
        parameter_names = ["log10_A"]

    elif kind == "smbhb_ptarcade_bhb_prior":
        mu, cov = bhb_priors["NG15"]
        sig = Ceffyl.signal(
            N_freqs=N_BINS,
            psd=powerlaw2_ceffyl,
            params=[parameter.Normal(mu=mu, sigma=cov, size=2)("gw_bhb")],
            name="",
        )
        parameter_names = ["gw_bhb_0", "gw_bhb_1"]
        chol = np.linalg.cholesky(cov)

        def prior_transform(u):
            uu = np.clip(np.asarray(u), 1e-12, 1.0 - 1e-12)
            return mu + chol @ ndtri(uu)

    elif kind in MODEL_FILES:
        model_path = MODEL_ROOT / MODEL_FILES[kind]
        mod = load_module(model_path, f"ng15_{kind}")
        sig = Ceffyl.signal(
            N_freqs=N_BINS,
            psd=aux.omega2cross(mod.spectrum, ceffyl=True),
            params=list(ParamDict(mod.parameters).values()),
            name="",
        )
        parameter_names = list(mod.parameters.keys())

    else:
        raise ValueError(f"unknown model kind: {kind}")

    cpta.add_signals([sig])
    if prior_transform is None:
        prior_transform = cpta.hypercube
    return cpta, parameter_names, prior_transform, model_path


def finite_loglike(cpta):
    def wrapped(theta):
        val = cpta.ln_likelihood(theta)
        if not np.isfinite(val):
            return -np.inf
        return float(val)

    return wrapped


def run_model(kind: str, cfg: dict):
    cpta, parameter_names, prior_transform, model_path = build_ceffyl_model(kind)
    t0 = time.time()
    sampler = NestedSampler(
        finite_loglike(cpta),
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
        "model": kind,
        "dim": int(cpta.ndim),
        "parameter_names": parameter_names,
        "ceffyl_param_names": list(cpta.param_names),
        "model_path": str(model_path) if model_path is not None else None,
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
    for cfg_label, models_out in config_results.items():
        if "smbhb_ptarcade_bhb_prior" not in models_out:
            continue
        base = models_out["smbhb_ptarcade_bhb_prior"]
        for kind in MODEL_FILES:
            if kind not in models_out:
                continue
            r = models_out[kind]
            lnB = float(r["lnZ"] - base["lnZ"])
            lnBerr = float(np.sqrt(r["lnZerr"] ** 2 + base["lnZerr"] ** 2))
            resid = float(lnB - PUBLISHED_LNB[kind])
            rows.append(
                {
                    "config": cfg_label,
                    "model": kind,
                    "lnB_vs_smbhb_ptarcade_bhb_prior": lnB,
                    "lnB_err": lnBerr,
                    "published_B": PUBLISHED_B[kind],
                    "published_lnB": PUBLISHED_LNB[kind],
                    "residual_vs_published_lnB": resid,
                }
            )

    aggregate = {}
    for kind in MODEL_FILES:
        vals = np.array(
            [r["lnB_vs_smbhb_ptarcade_bhb_prior"] for r in rows if r["model"] == kind],
            dtype=float,
        )
        if vals.size == 0:
            continue
        residuals = np.array(
            [r["residual_vs_published_lnB"] for r in rows if r["model"] == kind],
            dtype=float,
        )
        errs = np.array([r["lnB_err"] for r in rows if r["model"] == kind], dtype=float)
        aggregate[kind] = {
            "n_configs": int(vals.size),
            "lnB_mean": float(vals.mean()),
            "lnB_sample_std": float(vals.std(ddof=1)) if vals.size > 1 else 0.0,
            "lnB_min": float(vals.min()),
            "lnB_max": float(vals.max()),
            "mean_nested_lnB_err": float(errs.mean()),
            "published_lnB": PUBLISHED_LNB[kind],
            "residual_mean": float(residuals.mean()),
            "residual_sample_std": float(residuals.std(ddof=1)) if residuals.size > 1 else 0.0,
            "residual_min": float(residuals.min()),
            "residual_max": float(residuals.max()),
            "max_abs_residual": float(np.max(np.abs(residuals))),
        }

    return rows, aggregate


def write_output(path: Path, config_results: dict, complete: bool):
    rows, aggregate = summarize(config_results) if config_results else ([], {})
    out = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "complete": complete,
        "task": "PRL official-density robustness sweep",
        "source": {
            "ptarcade_package_version": ptarcade.__version__,
            "ptarcade_ceffyl_data_zenodo": "https://doi.org/10.5281/zenodo.10495907",
            "ng15_model_files_zenodo": "https://doi.org/10.5281/zenodo.8084351",
            "ceffyl_datadir": str(DATADIR),
            "model_root": str(MODEL_ROOT),
        },
        "method": {
            "mode": "dynesty evidence on PTArcade official ceffyl density grid",
            "n_bins": N_BINS,
            "configs": CONFIGS,
            "baseline": "smbhb_ptarcade_bhb_prior",
            "models": MODEL_FILES,
        },
        "config_results": config_results,
        "rows": rows,
        "aggregate": aggregate,
    }
    path.write_text(json.dumps(out, indent=2))


def main():
    print("=" * 72)
    print("PRL official-density robustness sweep")
    print("=" * 72)
    print(f"PTArcade package: {ptarcade.__version__}")
    print(f"Ceffyl datadir:   {DATADIR}")
    print(f"Models:           {', '.join(KINDS)}")
    config_results = {}

    for cfg in CONFIGS:
        label = str(cfg["label"])
        print(f"\n## config {label}: nlive={cfg['nlive']} dlogz={cfg['dlogz']} seed={cfg['seed']}")
        config_results[label] = {}
        for kind in KINDS:
            print(f"  running {kind} ...", flush=True)
            r = run_model(kind, cfg)
            config_results[label][kind] = r
            print(
                f"    lnZ={r['lnZ']:+.3f} +/- {r['lnZerr']:.3f}; "
                f"niter={r['niter']}; runtime={r['runtime_s']:.1f}s",
                flush=True,
            )
            write_output(PARTIAL, config_results, complete=False)

        rows, _ = summarize(config_results)
        recent = [r for r in rows if r["config"] == label]
        for row in recent:
            print(
                f"  {row['model']}: lnB={row['lnB_vs_smbhb_ptarcade_bhb_prior']:+.3f} "
                f"+/- {row['lnB_err']:.3f}; residual={row['residual_vs_published_lnB']:+.3f}",
                flush=True,
            )

    write_output(OUTFILE, config_results, complete=True)
    rows, aggregate = summarize(config_results)
    print("\n" + "=" * 72)
    print("Aggregate summary")
    print("=" * 72)
    for kind, agg in aggregate.items():
        print(
            f"{kind}: mean lnB={agg['lnB_mean']:+.3f}, "
            f"sample std={agg['lnB_sample_std']:.3f}, "
            f"mean residual={agg['residual_mean']:+.3f}, "
            f"max |residual|={agg['max_abs_residual']:.3f}"
        )
    print(f"saved: {OUTFILE}")


if __name__ == "__main__":
    main()
