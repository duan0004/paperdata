#!/usr/bin/env python3
"""PRL H9: cumulative-bin evidence driver on the official ceffyl density.

This is a cumulative low-frequency-bin scan.  For each N in N_BINS_LIST the
script evaluates the same stochastic spectral templates on the first N
Fourier bins of the official PTArcade ceffyl density.  It does not perform
timing-residual processing and does not touch pulsar par files.
"""

from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ceffyl import Ceffyl
from dynesty import NestedSampler
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
FIGDIR = ROOT / "results/T2_NG15yr/figures"
OUTFILE = OUTDIR / "prl_H9_bin_driver_analysis.json"
PARTIAL = OUTDIR / "prl_H9_bin_driver_analysis.partial.json"
OUT_MD = OUTDIR / "prl_H9_bin_driver_analysis.md"
OUT_FIG_PDF = FIGDIR / "prl_H9_bin_driver_figure.pdf"
OUT_FIG_PNG = FIGDIR / "prl_H9_bin_driver_figure.png"

N_BINS_LIST = [4, 6, 8, 10, 12, 14]
NLIVE = 500
DLOGZ = 0.1
SEED = 20260421

MODEL_ORDER = [
    "baseline_smbhb_ptarcade_bhb_prior",
    "sigw_gauss",
    "sigw_delta",
    "super",
    "ecc_supp_fixed_gamma",
    "broken_pl_fixed_gamma",
    "env_fixed_gamma",
]

MODEL_LABELS = {
    "sigw_gauss": "SIGW Gaussian",
    "sigw_delta": "SIGW delta",
    "super": "Cosmic superstrings",
    "ecc_supp_fixed_gamma": "Ecc.-inspired curvature",
    "broken_pl_fixed_gamma": "Broken-PL curvature",
    "env_fixed_gamma": "Environmental turnover",
}

MODEL_GROUP = {
    "sigw_gauss": "cosmology",
    "sigw_delta": "cosmology",
    "super": "cosmology",
    "ecc_supp_fixed_gamma": "astrophysical curvature",
    "broken_pl_fixed_gamma": "astrophysical curvature",
    "env_fixed_gamma": "astrophysical curvature",
}


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
    base = base_powerlaw_residual(f, Tspan, log10_A, gamma)
    fe = 10.0**log10_fe
    return base * (1.0 + (fe / f) ** 2.0) ** (-beta)


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


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


def baseline_bhb_prior(n_bins: int):
    cpta = Ceffyl.ceffyl(str(DATADIR))
    mu, cov = bhb_priors["NG15"]
    sig = Ceffyl.signal(
        N_freqs=n_bins,
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


def build_official_model(model_key: str, n_bins: int):
    path = MODEL_ROOT / f"{model_key}.py"
    mod = load_module(path, f"ng15_{model_key}_N{n_bins}")
    if not hasattr(mod, "spectrum"):
        raise ValueError(f"{model_key} does not define spectrum(f, ...)")
    cpta = Ceffyl.ceffyl(str(DATADIR))
    sig = Ceffyl.signal(
        N_freqs=n_bins,
        psd=aux.omega2cross(mod.spectrum, ceffyl=True),
        params=list(ParamDict(mod.parameters).values()),
        name="",
    )
    cpta.add_signals([sig])
    return cpta, list(mod.parameters.keys()), cpta.hypercube, path


def build_astro_model(model_key: str, n_bins: int):
    cpta = Ceffyl.ceffyl(str(DATADIR))
    if model_key == "env_fixed_gamma":
        params = [
            parameter.Uniform(-18, -11)("log10_A"),
            parameter.Uniform(-9.8, -7.0)("log10_fbend"),
            parameter.Uniform(0.5, 8.0)("kappa"),
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
    sig = Ceffyl.signal(N_freqs=n_bins, psd=psd, params=params, name="")
    cpta.add_signals([sig])
    return cpta, names, cpta.hypercube, None


def build_model(model_key: str, n_bins: int):
    if model_key == "baseline_smbhb_ptarcade_bhb_prior":
        return baseline_bhb_prior(n_bins)
    if model_key in {"sigw_gauss", "sigw_delta", "super"}:
        return build_official_model(model_key, n_bins)
    return build_astro_model(model_key, n_bins)


def run_model(model_key: str, n_bins: int):
    cpta, names, prior_transform, model_path = build_model(model_key, n_bins)
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
    best_idx = int(np.argmax(res["logl"]))
    return {
        "model": model_key,
        "N_bins": int(n_bins),
        "dim": int(cpta.ndim),
        "parameter_names": names,
        "ceffyl_param_names": list(cpta.param_names),
        "model_path": str(model_path) if model_path is not None else None,
        "lnZ": float(res["logz"][-1]),
        "lnZerr": float(res["logzerr"][-1]),
        "niter": int(res.niter),
        "runtime_s": float(time.time() - t0),
        "bestfit": res["samples"][best_idx].tolist(),
        "bestfit_loglike": float(res["logl"][best_idx]),
    }


def summarize(results_by_n: dict):
    rows = []
    for n_label, models in results_by_n.items():
        n_bins = int(n_label)
        base = models.get("baseline_smbhb_ptarcade_bhb_prior")
        if base is None:
            continue
        for key in MODEL_ORDER:
            if key == "baseline_smbhb_ptarcade_bhb_prior" or key not in models:
                continue
            r = models[key]
            rows.append(
                {
                    "N_bins": n_bins,
                    "model": key,
                    "group": MODEL_GROUP[key],
                    "lnB": float(r["lnZ"] - base["lnZ"]),
                    "lnB_err": float(np.sqrt(r["lnZerr"] ** 2 + base["lnZerr"] ** 2)),
                    "lnZ": r["lnZ"],
                    "lnZerr": r["lnZerr"],
                }
            )

    by_model = {}
    for key in MODEL_ORDER:
        if key == "baseline_smbhb_ptarcade_bhb_prior":
            continue
        vals = {row["N_bins"]: row for row in rows if row["model"] == key}
        if not vals:
            continue
        n_final = max(vals)
        final = vals[n_final]["lnB"]
        first_within_0p5 = None
        first_within_1p0 = None
        for n in sorted(vals):
            if first_within_1p0 is None and abs(vals[n]["lnB"] - final) <= 1.0:
                first_within_1p0 = n
            if first_within_0p5 is None and abs(vals[n]["lnB"] - final) <= 0.5:
                first_within_0p5 = n
        ordered_n = sorted(vals)
        if len(ordered_n) > 1:
            max_abs_step = max(
                abs(vals[n2]["lnB"] - vals[n1]["lnB"])
                for n1, n2 in zip(ordered_n[:-1], ordered_n[1:])
            )
        else:
            max_abs_step = 0.0
        by_model[key] = {
            "group": MODEL_GROUP[key],
            "lnB_final_N14": float(final),
            "first_N_within_1nat_of_N14": first_within_1p0,
            "first_N_within_0p5nat_of_N14": first_within_0p5,
            "delta_lnB_N4_to_N14": float(final - vals[4]["lnB"]) if 4 in vals else None,
            "delta_lnB_N6_to_N14": float(final - vals[6]["lnB"]) if 6 in vals else None,
            "delta_lnB_N8_to_N14": float(final - vals[8]["lnB"]) if 8 in vals else None,
            "max_abs_step": float(max_abs_step),
        }
    return rows, by_model


def write_json(path: Path, results_by_n: dict, complete: bool):
    rows, by_model = summarize(results_by_n)
    out = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "complete": complete,
        "task": "PRL H9 cumulative low-frequency-bin evidence driver",
        "source": {
            "ptarcade_package_version": ptarcade.__version__,
            "ceffyl_datadir": str(DATADIR),
            "model_root": str(MODEL_ROOT),
        },
        "method": {
            "n_bins_list": N_BINS_LIST,
            "nlive": NLIVE,
            "dlogz": DLOGZ,
            "seed": SEED,
            "baseline": "baseline_smbhb_ptarcade_bhb_prior",
            "scope": "cumulative first-N Fourier-bin scan on official PTArcade ceffyl density",
        },
        "dimension_checks": {
            "official_templates": "handled by ptarcade.models_utils.omega2cross with ceffyl=True",
            "astrophysical_curvature": "multipliers depend only on frequency ratios and are dimensionless",
            "low_frequency_limit": "curvature controls suppress or preserve residual power for non-negative shape parameters",
        },
        "results_by_n": results_by_n,
        "rows": rows,
        "summary_by_model": by_model,
    }
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")


def make_figure():
    data = json.loads(OUTFILE.read_text())
    rows = data["rows"]
    FIGDIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.0,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), constrained_layout=True)
    colors = {
        "sigw_gauss": "#26547c",
        "sigw_delta": "#3b73a3",
        "super": "#6c8ead",
        "ecc_supp_fixed_gamma": "#b36b00",
        "broken_pl_fixed_gamma": "#d08c13",
        "env_fixed_gamma": "#8a5a00",
    }
    markers = {
        "sigw_gauss": "o",
        "sigw_delta": "s",
        "super": "^",
        "ecc_supp_fixed_gamma": "D",
        "broken_pl_fixed_gamma": "P",
        "env_fixed_gamma": "v",
    }

    ax = axes[0]
    for key in MODEL_ORDER[1:]:
        key_rows = sorted([row for row in rows if row["model"] == key], key=lambda r: r["N_bins"])
        ax.errorbar(
            [row["N_bins"] for row in key_rows],
            [row["lnB"] for row in key_rows],
            yerr=[row["lnB_err"] for row in key_rows],
            marker=markers[key],
            color=colors[key],
            linewidth=1.2,
            markersize=3.5,
            capsize=2,
            label=MODEL_LABELS[key],
        )
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xlabel("first N low-frequency bins")
    ax.set_ylabel(r"$\ln B_N$ vs matched BHB prior")
    ax.set_title("A. Cumulative evidence")
    ax.grid(alpha=0.22, linewidth=0.6)
    ax.legend(loc="upper left", ncol=2, frameon=True, framealpha=0.88, fontsize=6.2)

    ax = axes[1]
    summary = data["summary_by_model"]
    ordered = sorted(summary.items(), key=lambda item: item[1]["lnB_final_N14"], reverse=True)
    y = np.arange(len(ordered))[::-1]
    vals = [item[1]["first_N_within_0p5nat_of_N14"] for item in ordered]
    bar_colors = ["#26547c" if item[1]["group"] == "cosmology" else "#b36b00" for item in ordered]
    ax.barh(y, vals, color=bar_colors, edgecolor="black", linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels([MODEL_LABELS[item[0]] for item in ordered])
    ax.set_xlim(0, 15)
    ax.set_xlabel(r"first N within 0.5 nat of $\ln B_{14}$")
    ax.set_title("B. Low-bin saturation")
    ax.grid(axis="x", alpha=0.22, linewidth=0.6)
    for ypos, val in zip(y, vals):
        ax.text(val + 0.25, ypos, str(val), va="center", ha="left", fontsize=7.5)

    fig.savefig(OUT_FIG_PDF)
    fig.savefig(OUT_FIG_PNG, dpi=250)


def write_md():
    data = json.loads(OUTFILE.read_text())
    lines = [
        "# PRL H9 Cumulative-Bin Evidence Driver",
        "",
        f"**Generated**: {data['generated']}",
        f"**JSON**: `{OUTFILE.relative_to(ROOT)}`",
        f"**Figure PDF**: `{OUT_FIG_PDF.relative_to(ROOT)}`",
        f"**Figure PNG**: `{OUT_FIG_PNG.relative_to(ROOT)}`",
        "",
        "## Summary",
        "",
        "| Model | group | lnB(N=14) | first N within 1 nat | first N within 0.5 nat | delta lnB N=6 to 14 |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for key, item in sorted(
        data["summary_by_model"].items(),
        key=lambda kv: kv[1]["lnB_final_N14"],
        reverse=True,
    ):
        lines.append(
            f"| `{key}` | {item['group']} | `{item['lnB_final_N14']:+.3f}` | "
            f"`{item['first_N_within_1nat_of_N14']}` | "
            f"`{item['first_N_within_0p5nat_of_N14']}` | "
            f"`{item['delta_lnB_N6_to_N14']:+.3f}` |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This cumulative-bin scan quantifies whether the matched-scale source",
            "ranking is already determined by the lowest Fourier bins.  It should be",
            "reported as a stochastic-template diagnostic, not as a time-domain",
            "source-population completeness claim.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main():
    print("PRL H9 cumulative-bin evidence driver")
    print(f"Ceffyl datadir: {DATADIR}")
    results_by_n = {}
    for n_bins in N_BINS_LIST:
        print(f"\n## N_bins={n_bins}", flush=True)
        results_by_n[str(n_bins)] = {}
        for key in MODEL_ORDER:
            print(f"running {key}", flush=True)
            r = run_model(key, n_bins)
            results_by_n[str(n_bins)][key] = r
            print(
                f"  lnZ={r['lnZ']:+.3f} +/- {r['lnZerr']:.3f}; "
                f"niter={r['niter']}; runtime={r['runtime_s']:.1f}s",
                flush=True,
            )
            write_json(PARTIAL, results_by_n, complete=False)
        rows, _ = summarize({str(n_bins): results_by_n[str(n_bins)]})
        for row in rows:
            print(
                f"  {row['model']}: lnB_N={row['lnB']:+.3f} +/- {row['lnB_err']:.3f}",
                flush=True,
            )
    write_json(OUTFILE, results_by_n, complete=True)
    make_figure()
    write_md()
    print(f"saved: {OUTFILE}")
    print(f"saved: {OUT_MD}")
    print(f"saved: {OUT_FIG_PDF}")
    print(f"saved: {OUT_FIG_PNG}")


if __name__ == "__main__":
    main()
