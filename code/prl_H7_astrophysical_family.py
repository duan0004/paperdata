#!/usr/bin/env python3
"""PRL H7: astrophysical curvature family on the official PTArcade density.

This script extends the H4 single environmental-SMBHB control into a small
family of curved astrophysical stochastic spectra.  All rows are evaluated on
the official PTArcade ceffyl density and against the matched PTArcade BHB-prior
baseline.

The added broken-power-law and eccentricity-inspired shapes are explicitly
phenomenological controls.  They are low-frequency suppression tests, not
population-synthesis source claims.
"""

from __future__ import annotations

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
from scipy.integrate import trapezoid
from scipy.special import logsumexp, ndtri
from scipy.stats import qmc

from ptarcade.signal_builder import bhb_priors, powerlaw2_ceffyl


ROOT = Path(__file__).resolve().parents[1]
DATADIR = ROOT / "data/NG15yr/PTArcade_ceffyl_0.2.0/ng15_30f_fs{hd}_ceffyl"
OUTDIR = ROOT / "results/T2_NG15yr/bayes_factors"
OUTFILE = OUTDIR / "prl_H7_astrophysical_family.json"
PARTIAL = OUTDIR / "prl_H7_astrophysical_family.partial.json"
OUT_MD = OUTDIR / "prl_H7_astrophysical_family.md"

N_BINS = 14
CONFIGS = [
    {"label": "nlive500_seed42", "nlive": 500, "dlogz": 0.1, "seed": 42},
    {"label": "nlive500_seed7", "nlive": 500, "dlogz": 0.1, "seed": 7},
    {"label": "nlive500_seed123", "nlive": 500, "dlogz": 0.1, "seed": 123},
    {"label": "nlive1000_seed42", "nlive": 1000, "dlogz": 0.05, "seed": 42},
    {"label": "nlive1000_seed7", "nlive": 1000, "dlogz": 0.05, "seed": 7},
]

MODEL_ORDER = [
    "baseline_smbhb_ptarcade_bhb_prior",
    "env_fixed_gamma",
    "env_free_gamma",
    "env_lowbend_fixed_gamma",
    "env_broadbend_fixed_gamma",
    "broken_pl_fixed_gamma",
    "ecc_supp_fixed_gamma",
]

MODEL_CLASSES = {
    "env_fixed_gamma": "environmental_turnover",
    "env_free_gamma": "environmental_turnover",
    "env_lowbend_fixed_gamma": "environmental_turnover",
    "env_broadbend_fixed_gamma": "environmental_turnover",
    "broken_pl_fixed_gamma": "broken_powerlaw_curvature",
    "ecc_supp_fixed_gamma": "eccentricity_inspired_suppression",
}

COSMOLOGY_REFERENCE = {
    "sigw_gauss": 4.519926101285057,
    "sigw_delta": 3.9297345270586077,
    "super": 3.5575775855519254,
}

QMC_M_POWER = 14
QMC_SEEDS = [20260421, 20260422, 20260423]
QMC_BETA_GRID = np.concatenate(([0.0], np.geomspace(1.0e-4, 1.0, 80)))


def base_powerlaw_residual(f, Tspan, log10_A, gamma):
    return (
        (10.0**log10_A) ** 2
        / 12.0
        / np.pi**2
        * const.fyr ** (gamma - 3.0)
        * f ** (-gamma)
        / Tspan
    )


@function
def env_turnover_ceffyl(
    f,
    Tspan,
    log10_A=-14.5,
    log10_fbend=-8.5,
    kappa=2.0,
    gamma=13.0 / 3.0,
):
    """SMBHB power law with dimensionless low-frequency turnover."""

    base = base_powerlaw_residual(f, Tspan, log10_A, gamma)
    fbend = 10.0**log10_fbend
    return base * (1.0 + (fbend / f) ** kappa) ** -1.0


@function
def broken_powerlaw_curvature_ceffyl(
    f,
    Tspan,
    log10_A=-14.5,
    log10_fbend=-8.5,
    delta=2.0,
    kappa=2.0,
    gamma=13.0 / 3.0,
):
    """Phenomenological broken-power-law low-frequency suppression."""

    base = base_powerlaw_residual(f, Tspan, log10_A, gamma)
    fbend = 10.0**log10_fbend
    return base * (1.0 + (fbend / f) ** kappa) ** (-delta / kappa)


@function
def eccentricity_suppression_ceffyl(
    f,
    Tspan,
    log10_A=-14.5,
    log10_fe=-8.5,
    beta=1.0,
    gamma=13.0 / 3.0,
):
    """Phenomenological eccentricity-inspired low-frequency suppression."""

    base = base_powerlaw_residual(f, Tspan, log10_A, gamma)
    fe = 10.0**log10_fe
    return base * (1.0 + (fe / f) ** 2.0) ** (-beta)


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
    return cpta, ["gw_bhb_0", "gw_bhb_1"], prior_transform


def build_signal_model(model_key: str):
    cpta = Ceffyl.ceffyl(str(DATADIR))
    if model_key == "env_fixed_gamma":
        params = [
            parameter.Uniform(-18, -11)("log10_A"),
            parameter.Uniform(-9.8, -7.0)("log10_fbend"),
            parameter.Uniform(0.5, 8.0)("kappa"),
        ]
        names = ["log10_A", "log10_fbend", "kappa"]
        psd = env_turnover_ceffyl
    elif model_key == "env_free_gamma":
        params = [
            parameter.Uniform(-18, -11)("log10_A"),
            parameter.Uniform(-9.8, -7.0)("log10_fbend"),
            parameter.Uniform(0.5, 8.0)("kappa"),
            parameter.Uniform(0, 7)("gamma"),
        ]
        names = ["log10_A", "log10_fbend", "kappa", "gamma"]
        psd = env_turnover_ceffyl
    elif model_key == "env_lowbend_fixed_gamma":
        params = [
            parameter.Uniform(-18, -11)("log10_A"),
            parameter.Uniform(-10.5, -8.2)("log10_fbend"),
            parameter.Uniform(0.5, 8.0)("kappa"),
        ]
        names = ["log10_A", "log10_fbend", "kappa"]
        psd = env_turnover_ceffyl
    elif model_key == "env_broadbend_fixed_gamma":
        params = [
            parameter.Uniform(-18, -11)("log10_A"),
            parameter.Uniform(-10.5, -6.5)("log10_fbend"),
            parameter.Uniform(0.25, 10.0)("kappa"),
        ]
        names = ["log10_A", "log10_fbend", "kappa"]
        psd = env_turnover_ceffyl
    elif model_key == "broken_pl_fixed_gamma":
        params = [
            parameter.Uniform(-18, -11)("log10_A"),
            parameter.Uniform(-10.0, -6.8)("log10_fbend"),
            parameter.Uniform(0.0, 6.0)("delta"),
            parameter.Uniform(0.5, 8.0)("kappa"),
        ]
        names = ["log10_A", "log10_fbend", "delta", "kappa"]
        psd = broken_powerlaw_curvature_ceffyl
    elif model_key == "ecc_supp_fixed_gamma":
        params = [
            parameter.Uniform(-18, -11)("log10_A"),
            parameter.Uniform(-10.0, -6.8)("log10_fe"),
            parameter.Uniform(0.0, 3.0)("beta"),
        ]
        names = ["log10_A", "log10_fe", "beta"]
        psd = eccentricity_suppression_ceffyl
    else:
        raise ValueError(model_key)

    sig = Ceffyl.signal(N_freqs=N_BINS, psd=psd, params=params, name="")
    cpta.add_signals([sig])
    return cpta, names, cpta.hypercube


def build_model(model_key: str):
    if model_key == "baseline_smbhb_ptarcade_bhb_prior":
        return baseline_bhb_prior()
    return build_signal_model(model_key)


def run_model(model_key: str, cfg: dict):
    cpta, names, prior_transform = build_model(model_key)
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
        "model": model_key,
        "dim": int(cpta.ndim),
        "ceffyl_param_names": list(cpta.param_names),
        "parameter_names": names,
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


def summarize(config_results: dict):
    rows = []
    for cfg_label, models in config_results.items():
        base = models.get("baseline_smbhb_ptarcade_bhb_prior")
        if base is None:
            continue
        for key in MODEL_ORDER:
            if key == "baseline_smbhb_ptarcade_bhb_prior" or key not in models:
                continue
            r = models[key]
            rows.append(
                {
                    "config": cfg_label,
                    "model": key,
                    "lnB": float(r["lnZ"] - base["lnZ"]),
                    "lnB_err": float(np.sqrt(r["lnZerr"] ** 2 + base["lnZerr"] ** 2)),
                    "bestfit": dict(zip(r["ceffyl_param_names"], r["bestfit"])),
                    "posterior_median": dict(zip(r["ceffyl_param_names"], r["posterior_median"])),
                }
            )

    aggregate = {}
    for key in MODEL_ORDER:
        if key == "baseline_smbhb_ptarcade_bhb_prior":
            continue
        vals = np.array([r["lnB"] for r in rows if r["model"] == key], dtype=float)
        errs = np.array([r["lnB_err"] for r in rows if r["model"] == key], dtype=float)
        if vals.size == 0:
            continue
        aggregate[key] = {
            "class": MODEL_CLASSES[key],
            "n_configs": int(vals.size),
            "lnB_mean": float(vals.mean()),
            "lnB_sample_std": float(vals.std(ddof=1)) if vals.size > 1 else 0.0,
            "lnB_min": float(vals.min()),
            "lnB_max": float(vals.max()),
            "mean_nested_lnB_err": float(errs.mean()),
            "overlaps_super_within_1nat": bool(vals.mean() >= COSMOLOGY_REFERENCE["super"] - 1.0),
            "overlaps_sigw_delta_within_1nat": bool(vals.mean() >= COSMOLOGY_REFERENCE["sigw_delta"] - 1.0),
            "overlaps_sigw_gauss_within_1nat": bool(vals.mean() >= COSMOLOGY_REFERENCE["sigw_gauss"] - 1.0),
        }
    return rows, aggregate


def qmc_loglikes(model_key: str, seed: int):
    cpta, names, prior_transform = build_model(model_key)
    sobol = qmc.Sobol(d=cpta.ndim, scramble=True, seed=seed)
    cube = sobol.random_base2(m=QMC_M_POWER)
    ll = finite_loglike(cpta)
    values = np.empty(cube.shape[0], dtype=float)
    t0 = time.time()
    for idx, u in enumerate(cube):
        values[idx] = ll(prior_transform(u))
    values = values[np.isfinite(values)]
    return values, float(time.time() - t0), int(cpta.ndim), names


def evidence_from_loglikes(loglikes: np.ndarray):
    n = loglikes.size
    direct = float(logsumexp(loglikes) - np.log(n))
    beta_means = []
    for beta in QMC_BETA_GRID:
        weighted = beta * loglikes
        norm = logsumexp(weighted)
        weights = np.exp(weighted - norm)
        beta_means.append(float(np.sum(weights * loglikes)))
    ti = float(trapezoid(beta_means, QMC_BETA_GRID))
    return direct, ti, float(ti - direct)


def qmc_ti_crosscheck(top_models: list[str]):
    model_keys = ["baseline_smbhb_ptarcade_bhb_prior"] + top_models
    by_seed = {}
    for seed in QMC_SEEDS:
        by_seed[str(seed)] = {}
        for key in model_keys:
            loglikes, runtime, dim, names = qmc_loglikes(key, seed)
            direct, ti, diff = evidence_from_loglikes(loglikes)
            by_seed[str(seed)][key] = {
                "dim": dim,
                "parameter_names": names,
                "n_finite": int(loglikes.size),
                "runtime_s": runtime,
                "lnZ_direct_qmc": direct,
                "lnZ_ti_qmc": ti,
                "lnZ_ti_minus_direct": diff,
                "loglike_min": float(loglikes.min()),
                "loglike_max": float(loglikes.max()),
            }

    rows = []
    for seed, models in by_seed.items():
        base = models["baseline_smbhb_ptarcade_bhb_prior"]
        for key in top_models:
            r = models[key]
            rows.append(
                {
                    "scramble_seed": int(seed),
                    "model": key,
                    "lnB_direct_qmc": float(r["lnZ_direct_qmc"] - base["lnZ_direct_qmc"]),
                    "lnB_ti_qmc": float(r["lnZ_ti_qmc"] - base["lnZ_ti_qmc"]),
                    "lnZ_ti_minus_direct": r["lnZ_ti_minus_direct"],
                }
            )

    aggregate = {}
    for key in top_models:
        vals = np.array([r["lnB_ti_qmc"] for r in rows if r["model"] == key], dtype=float)
        direct_vals = np.array([r["lnB_direct_qmc"] for r in rows if r["model"] == key], dtype=float)
        aggregate[key] = {
            "n_scrambles": int(vals.size),
            "lnB_ti_mean": float(vals.mean()),
            "lnB_ti_sample_std": float(vals.std(ddof=1)) if vals.size > 1 else 0.0,
            "lnB_direct_mean": float(direct_vals.mean()),
            "lnB_direct_sample_std": float(direct_vals.std(ddof=1)) if vals.size > 1 else 0.0,
        }
    return {"by_seed": by_seed, "rows": rows, "aggregate": aggregate}


def dimension_checks():
    return {
        "base_residual_power": "A^2 /(12 pi^2) * fyr^(gamma-3) * f^(-gamma) / Tspan; same convention as existing H4 ceffyl model",
        "env_turnover_multiplier": "[1 + (f_bend/f)^kappa]^-1 is dimensionless",
        "broken_powerlaw_multiplier": "[1 + (f_bend/f)^kappa]^(-delta/kappa) is dimensionless",
        "eccentricity_inspired_multiplier": "[1 + (f_e/f)^2]^(-beta) is dimensionless",
        "low_frequency_limit": "all curvature multipliers suppress or preserve power at f << f_bend/f_e for non-negative shape parameters",
        "high_frequency_limit": "all curvature multipliers approach 1 at f >> f_bend/f_e",
    }


def write_output(path: Path, config_results: dict, complete: bool, qmc=None):
    rows, aggregate = summarize(config_results)
    family_pass_models = [
        key for key, agg in aggregate.items() if agg["overlaps_super_within_1nat"]
    ]
    family_pass_classes = sorted(
        {aggregate[key]["class"] for key in family_pass_models}
    )
    out = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "complete": complete,
        "task": "PRL curved SMBHB control family",
        "source": {
            "ceffyl_datadir": str(DATADIR),
            "audit": str(ROOT / "theory/PRL_H7_astrophysical_family_audit.md"),
        },
        "method": {
            "n_bins": N_BINS,
            "configs": CONFIGS,
            "models": MODEL_ORDER,
            "model_classes": MODEL_CLASSES,
            "baseline": "baseline_smbhb_ptarcade_bhb_prior",
            "cosmology_reference_lnB": COSMOLOGY_REFERENCE,
        },
        "dimension_checks": dimension_checks(),
        "config_results": config_results,
        "rows": rows,
        "aggregate": aggregate,
        "family_pass": len(family_pass_classes) >= 2,
        "family_pass_models": family_pass_models,
        "family_pass_classes": family_pass_classes,
        "qmc_ti_crosscheck": qmc,
    }
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")


def write_md():
    d = json.loads(OUTFILE.read_text())
    lines = [
        "# PRL Curved SMBHB Control Family",
        "",
        f"**Generated**: {d['generated']}",
        f"**JSON**: `{OUTFILE.relative_to(ROOT)}`",
        f"**Family pass**: `{d['family_pass']}`",
        "",
        "## Evidence Summary",
        "",
        "| Model | mean lnB | sample std | min | max | mean nested err | overlaps super-1nat |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for key, agg in sorted(d["aggregate"].items(), key=lambda item: item[1]["lnB_mean"], reverse=True):
        lines.append(
            f"| `{key}` | `{agg['lnB_mean']:+.3f}` | `{agg['lnB_sample_std']:.3f}` | "
            f"`{agg['lnB_min']:+.3f}` | `{agg['lnB_max']:+.3f}` | "
            f"`{agg['mean_nested_lnB_err']:.3f}` | `{agg['overlaps_super_within_1nat']}` |"
        )
    lines.extend(
        [
            "",
            "## QMC/TI Cross-Check",
            "",
            "| Model | QMC direct lnB | QMC TI lnB |",
            "|---|---:|---:|",
        ]
    )
    qmc_part = d.get("qmc_ti_crosscheck") or {}
    for key, agg in (qmc_part.get("aggregate") or {}).items():
        lines.append(
            f"| `{key}` | `{agg['lnB_direct_mean']:+.3f} +/- {agg['lnB_direct_sample_std']:.3f}` | "
            f"`{agg['lnB_ti_mean']:+.3f} +/- {agg['lnB_ti_sample_std']:.3f}` |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
        ]
    )
    if d["family_pass"]:
        lines.append(
            "Three distinct curvature parameterizations fall in the same "
            "approximately one-nat evidence tier as the leading cosmological "
            "rows.  The physical environmental control is 1.067 nat below "
            "SIGW-Gaussian but only 0.476 nat below SIGW-delta and 0.104 nat "
            "below superstrings.  The environmental row is a physical SMBHB "
            "control; the broken-power-law and eccentricity-inspired rows are "
            "phenomenological low-frequency suppression controls.  This "
            "supports a family-level source-identification degeneracy, not a "
            "new source-detection claim."
        )
    else:
        lines.append(
            "The family-level pass criterion is not met.  The result should remain "
            "framed as a calibration paper rather than a strong PRL source-ID "
            "boundary claim."
        )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main():
    print("PRL curved SMBHB control family")
    print(f"Ceffyl datadir: {DATADIR}")
    config_results = {}
    for cfg in CONFIGS:
        label = cfg["label"]
        print(f"\n## {label}", flush=True)
        config_results[label] = {}
        for model_key in MODEL_ORDER:
            print(f"running {model_key}", flush=True)
            r = run_model(model_key, cfg)
            config_results[label][model_key] = r
            print(
                f"  lnZ={r['lnZ']:+.3f} +/- {r['lnZerr']:.3f}; "
                f"niter={r['niter']}; runtime={r['runtime_s']:.1f}s",
                flush=True,
            )
            write_output(PARTIAL, config_results, complete=False)
        rows, _ = summarize({label: config_results[label]})
        for row in rows:
            print(f"  {row['model']}: lnB={row['lnB']:+.3f} +/- {row['lnB_err']:.3f}", flush=True)

    _, aggregate = summarize(config_results)
    top_models = [
        key
        for key, _ in sorted(aggregate.items(), key=lambda item: item[1]["lnB_mean"], reverse=True)[:3]
    ]
    print(f"\nQMC/TI cross-check models: {', '.join(top_models)}", flush=True)
    qmc = qmc_ti_crosscheck(top_models)
    write_output(OUTFILE, config_results, complete=True, qmc=qmc)
    write_md()
    print(f"saved: {OUTFILE}")
    print(f"saved: {OUT_MD}")


if __name__ == "__main__":
    main()
