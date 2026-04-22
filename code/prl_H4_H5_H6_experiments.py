#!/usr/bin/env python3
"""
PRL H4/H5/H6 experiments.

H4: astrophysical curved-SMBHB competitor against the same official density.
H5: official-density evidence extension to all stochastic-background PTArcade
    model files available in the downloaded NG15 model package.
H6: SIGW-delta prior-boundary sensitivity for the peak-frequency upper bound.

Outputs:
  results/T2_NG15yr/bayes_factors/prl_H4_H5_H6_experiments.json
  results/T2_NG15yr/bayes_factors/prl_H4_H5_H6_experiments.md

Important scope:
  The script only evaluates PTArcade models that expose spectrum(f, ...).
  Models that expose signal(toas, ...) such as PBH/ULDM timing-residual
  signals are deliberately skipped because the released free-spectrum ceffyl
  density is not the correct likelihood object for them.
"""

from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from ceffyl import Ceffyl
from dynesty import NestedSampler
from dynesty.utils import resample_equal
import enterprise.constants as const
from enterprise.signals import parameter
from enterprise.signals.parameter import function
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
OUTFILE = OUTDIR / "prl_H4_H5_H6_experiments.json"
PARTIAL = OUTDIR / "prl_H4_H5_H6_experiments.partial.json"
OUT_MD = OUTDIR / "prl_H4_H5_H6_experiments.md"

N_BINS = 14
NLIVE = 500
DLOGZ = 0.1
SEED = 20260421

OFFICIAL_SPECTRUM_MODELS = [
    "dw_ds",
    "dw_sm",
    "igw",
    "meta_l",
    "meta_ls",
    "pt_bubble",
    "pt_sound",
    "sigw_box",
    "sigw_delta",
    "sigw_gauss",
    "stable_c",
    "stable_k",
    "stable_m",
    "stable_n",
    "super",
]

SKIPPED_SIGNAL_MODELS = [
    "pbh_dynamic",
    "pbh_static",
    "uldm_c_cor",
    "uldm_c_unc",
    "uldm_e",
    "uldm_p_cor",
    "uldm_p_unc",
    "uldm_vecBL_cor",
    "uldm_vecBL_unc",
    "uldm_vecB_cor",
    "uldm_vecB_unc",
]

PUBLISHED_LNB = {
    "sigw_gauss": float(np.log(57.0)),
    "sigw_delta": float(np.log(44.0)),
    "super": float(np.log(46.0)),
    "pt_bubble": float(np.log(18.0)),
    "pt_sound": float(np.log(3.7)),
}


@function
def smbhb_env_turnover_ceffyl(
    f,
    Tspan,
    log10_A=-14.5,
    log10_fbend=-8.5,
    kappa=2.0,
    gamma=13.0 / 3.0,
):
    """Phenomenological low-frequency turnover SMBHB spectrum.

    This is a curved-SMBHB control model, not a final population
    synthesis model.  The high-frequency limit is the standard PTA power law.
    At f << f_bend, the residual power is suppressed by
    [1 + (f_bend/f)^kappa]^{-1}.
    """

    base = (
        (10.0**log10_A) ** 2
        / 12.0
        / np.pi**2
        * const.fyr ** (gamma - 3.0)
        * f ** (-gamma)
        / Tspan
    )
    fbend = 10.0**log10_fbend
    turnover = (1.0 + (fbend / f) ** kappa) ** -1.0
    return base * turnover


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def baseline_bhb_prior():
    cpta = Ceffyl.ceffyl(str(DATADIR))
    mu, cov = bhb_priors["NG15"]
    sig = Ceffyl.signal(
        N_freqs=N_BINS,
        psd=powerlaw2_ceffyl,
        params=[parameter.Normal(mu=mu, sigma=cov, size=2)("gw_bhb")],
        name="",
    )
    chol = np.linalg.cholesky(cov)

    def prior_transform(u):
        uu = np.clip(np.asarray(u), 1e-12, 1.0 - 1e-12)
        return mu + chol @ ndtri(uu)

    cpta.add_signals([sig])
    return cpta, ["gw_bhb_0", "gw_bhb_1"], prior_transform, None


def build_env_model(kind: str):
    cpta = Ceffyl.ceffyl(str(DATADIR))
    if kind == "smbhb_env_fixed_gamma":
        params = [
            parameter.Uniform(-18, -11)("log10_A"),
            parameter.Uniform(-9.8, -7.0)("log10_fbend"),
            parameter.Uniform(0.5, 8.0)("kappa"),
        ]
        names = ["log10_A", "log10_fbend", "kappa"]
    elif kind == "smbhb_env_free_gamma":
        params = [
            parameter.Uniform(-18, -11)("log10_A"),
            parameter.Uniform(-9.8, -7.0)("log10_fbend"),
            parameter.Uniform(0.5, 8.0)("kappa"),
            parameter.Uniform(0, 7)("gamma"),
        ]
        names = ["log10_A", "log10_fbend", "kappa", "gamma"]
    else:
        raise ValueError(kind)
    sig = Ceffyl.signal(
        N_freqs=N_BINS,
        psd=smbhb_env_turnover_ceffyl,
        params=params,
        name="",
    )
    cpta.add_signals([sig])
    return cpta, names, cpta.hypercube, None


def build_official_spectrum_model(model_name: str):
    path = MODEL_ROOT / f"{model_name}.py"
    mod = load_module(path, f"ng15_{model_name}")
    if not hasattr(mod, "spectrum"):
        raise ValueError(f"{model_name} does not define spectrum(f, ...)")
    cpta = Ceffyl.ceffyl(str(DATADIR))
    sig = Ceffyl.signal(
        N_freqs=N_BINS,
        psd=aux.omega2cross(mod.spectrum, ceffyl=True),
        params=list(ParamDict(mod.parameters).values()),
        name="",
    )
    cpta.add_signals([sig])
    return cpta, list(mod.parameters.keys()), cpta.hypercube, path


def build_sigw_delta_prior_model(label: str, f_upper: float):
    path = MODEL_ROOT / "sigw_delta.py"
    mod = load_module(path, f"ng15_sigw_delta_{label}")
    cpta = Ceffyl.ceffyl(str(DATADIR))
    params = [
        parameter.Uniform(-11, f_upper)("log10_f_peak"),
        parameter.Uniform(-3, 1)("log10_A"),
    ]
    sig = Ceffyl.signal(
        N_freqs=N_BINS,
        psd=aux.omega2cross(mod.spectrum, ceffyl=True),
        params=params,
        name="",
    )
    cpta.add_signals([sig])
    return cpta, ["log10_f_peak", "log10_A"], cpta.hypercube, path


def finite_loglike(cpta):
    def wrapped(theta):
        try:
            val = cpta.ln_likelihood(theta)
        except Exception:
            return -np.inf
        if not np.isfinite(val):
            return -np.inf
        return float(val)

    return wrapped


def run_model(label: str, builder):
    cpta, parameter_names, prior_transform, model_path = builder()
    print(f"\n==== {label} (dim={cpta.ndim}) ====", flush=True)
    if model_path is not None:
        print(f"model file: {model_path}", flush=True)
    t0 = time.time()
    sampler = NestedSampler(
        finite_loglike(cpta),
        prior_transform,
        cpta.ndim,
        nlive=NLIVE,
        rstate=np.random.default_rng(SEED),
        bound="multi",
        sample="rwalk",
    )
    sampler.run_nested(dlogz=DLOGZ, print_progress=False)
    res = sampler.results
    weights = np.exp(res["logwt"] - res["logz"][-1])
    equal = resample_equal(res["samples"], weights)
    best_idx = int(np.argmax(res["logl"]))
    out = {
        "label": label,
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
        "posterior_median": np.median(equal, axis=0).tolist(),
        "posterior_p16": np.percentile(equal, 16, axis=0).tolist(),
        "posterior_p84": np.percentile(equal, 84, axis=0).tolist(),
    }
    print(
        f"lnZ={out['lnZ']:+.3f} +/- {out['lnZerr']:.3f}; "
        f"niter={out['niter']}; runtime={out['runtime_s']:.1f}s",
        flush=True,
    )
    return out


def attach_lnb(results: dict):
    base = results["baseline_smbhb_ptarcade_bhb_prior"]
    for label, r in results.items():
        r["lnB_vs_smbhb_ptarcade_bhb_prior"] = float(r["lnZ"] - base["lnZ"])
        r["lnB_err_vs_smbhb_ptarcade_bhb_prior"] = float(
            np.sqrt(r["lnZerr"] ** 2 + base["lnZerr"] ** 2)
        )
        if label.startswith("official_"):
            model = label.removeprefix("official_")
            if model in PUBLISHED_LNB:
                r["published_lnB"] = PUBLISHED_LNB[model]
                r["residual_vs_published_lnB"] = float(
                    r["lnB_vs_smbhb_ptarcade_bhb_prior"] - PUBLISHED_LNB[model]
                )


def summarize(results: dict):
    h4 = {
        k: results[k]
        for k in ["smbhb_env_fixed_gamma", "smbhb_env_free_gamma"]
        if k in results
    }
    official_rows = []
    for k, r in results.items():
        if k.startswith("official_"):
            model = k.removeprefix("official_")
            official_rows.append(
                {
                    "model": model,
                    "dim": r["dim"],
                    "lnB": r["lnB_vs_smbhb_ptarcade_bhb_prior"],
                    "lnB_err": r["lnB_err_vs_smbhb_ptarcade_bhb_prior"],
                    "published_lnB": r.get("published_lnB"),
                    "residual_vs_published_lnB": r.get("residual_vs_published_lnB"),
                    "bestfit": dict(zip(r["ceffyl_param_names"], r["bestfit"])),
                }
            )
    official_rows.sort(key=lambda x: x["lnB"], reverse=True)

    h6_rows = []
    original = results.get("sigw_delta_prior_fmax_-5")
    for k, r in results.items():
        if k.startswith("sigw_delta_prior_"):
            h6_rows.append(
                {
                    "label": k,
                    "lnB": r["lnB_vs_smbhb_ptarcade_bhb_prior"],
                    "lnB_err": r["lnB_err_vs_smbhb_ptarcade_bhb_prior"],
                    "delta_lnB_vs_original": None
                    if original is None
                    else float(
                        r["lnB_vs_smbhb_ptarcade_bhb_prior"]
                        - original["lnB_vs_smbhb_ptarcade_bhb_prior"]
                    ),
                    "bestfit": dict(zip(r["ceffyl_param_names"], r["bestfit"])),
                }
            )
    h6_rows.sort(key=lambda x: x["label"])
    return {"H4_env_smbhb": h4, "H5_official_rank": official_rows, "H6_sigw_delta_prior": h6_rows}


def write_json(path: Path, results: dict, complete: bool):
    attach_lnb(results)
    out = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "complete": complete,
        "task": "PRL H4/H5/H6 experiments",
        "source": {
            "ptarcade_package_version": ptarcade.__version__,
            "ptarcade_ceffyl_data_zenodo": "https://doi.org/10.5281/zenodo.10495907",
            "ng15_model_files_zenodo": "https://doi.org/10.5281/zenodo.8084351",
            "ceffyl_datadir": str(DATADIR),
            "model_root": str(MODEL_ROOT),
        },
        "method": {
            "n_bins": N_BINS,
            "nlive": NLIVE,
            "dlogz": DLOGZ,
            "seed": SEED,
            "baseline": "smbhb_ptarcade_bhb_prior",
            "skipped_signal_models": SKIPPED_SIGNAL_MODELS,
        },
        "results": results,
        "summary": summarize(results),
    }
    path.write_text(json.dumps(out, indent=2))


def write_markdown(json_path: Path):
    d = json.loads(json_path.read_text())
    summary = d["summary"]
    lines = [
        "# PRL H4/H5/H6 Experiment Results",
        "",
        f"**Date**: {d['generated']}",
        f"**JSON**: `{json_path.relative_to(ROOT)}`",
        f"**Baseline**: `{d['method']['baseline']}`",
        "",
        "## H4: Astrophysical Curved-SMBHB Competitor",
        "",
        "| Model | lnB vs BHB prior | lnB err | Best-fit parameters |",
        "|---|---:|---:|---|",
    ]
    for label, r in summary["H4_env_smbhb"].items():
        best = dict(zip(r["ceffyl_param_names"], r["bestfit"]))
        lines.append(
            f"| `{label}` | `{r['lnB_vs_smbhb_ptarcade_bhb_prior']:+.3f}` | "
            f"`{r['lnB_err_vs_smbhb_ptarcade_bhb_prior']:.3f}` | `{best}` |"
        )

    lines += [
        "",
        "## H5: Official-Density Stochastic-Background Model Ranking",
        "",
        "Only PTArcade model files defining `spectrum(f, ...)` were evaluated.",
        "PBH and ULDM `signal(toas, ...)` models were skipped by design.",
        "",
        "| Rank | Model | dim | lnB vs BHB prior | lnB err | published lnB | residual |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ]
    for i, row in enumerate(summary["H5_official_rank"], 1):
        pub = "" if row["published_lnB"] is None else f"`{row['published_lnB']:+.3f}`"
        resid = (
            ""
            if row["residual_vs_published_lnB"] is None
            else f"`{row['residual_vs_published_lnB']:+.3f}`"
        )
        lines.append(
            f"| {i} | `{row['model']}` | {row['dim']} | `{row['lnB']:+.3f}` | "
            f"`{row['lnB_err']:.3f}` | {pub} | {resid} |"
        )

    lines += [
        "",
        "## H6: SIGW-Delta Prior-Boundary Sensitivity",
        "",
        "| Prior label | lnB vs BHB prior | lnB err | delta lnB vs original | Best-fit parameters |",
        "|---|---:|---:|---:|---|",
    ]
    for row in summary["H6_sigw_delta_prior"]:
        delta = "" if row["delta_lnB_vs_original"] is None else f"`{row['delta_lnB_vs_original']:+.3f}`"
        lines.append(
            f"| `{row['label']}` | `{row['lnB']:+.3f}` | `{row['lnB_err']:.3f}` | "
            f"{delta} | `{row['bestfit']}` |"
        )

    lines += [
        "",
        "## Interpretation Placeholder",
        "",
        "Use this table to decide whether the PRL claim remains merely a provenance",
        "audit or becomes a broader source-identification systematic: if the curved",
        "SMBHB model is competitive with the leading cosmological templates, the PRL",
        "story should foreground curved-SMBHB source-identification degeneracy.",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


def main():
    print("=" * 72)
    print("PRL H4/H5/H6 experiments")
    print("=" * 72)
    print(f"PTArcade package: {ptarcade.__version__}")
    print(f"Ceffyl datadir: {DATADIR}")
    results = {}

    jobs = [
        ("baseline_smbhb_ptarcade_bhb_prior", baseline_bhb_prior),
        ("smbhb_env_fixed_gamma", lambda: build_env_model("smbhb_env_fixed_gamma")),
        ("smbhb_env_free_gamma", lambda: build_env_model("smbhb_env_free_gamma")),
    ]

    for model in OFFICIAL_SPECTRUM_MODELS:
        jobs.append((f"official_{model}", lambda model=model: build_official_spectrum_model(model)))

    for upper in [-5, -4, -3]:
        jobs.append(
            (
                f"sigw_delta_prior_fmax_{upper}",
                lambda upper=upper: build_sigw_delta_prior_model(f"fmax_{upper}", upper),
            )
        )

    for label, builder in jobs:
        results[label] = run_model(label, builder)
        write_json(PARTIAL, results, complete=False)

    write_json(OUTFILE, results, complete=True)
    write_markdown(OUTFILE)
    print("\n" + "=" * 72)
    print("Saved")
    print("=" * 72)
    print(OUTFILE)
    print(OUT_MD)


if __name__ == "__main__":
    main()
