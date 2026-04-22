#!/usr/bin/env python3
"""Independent evidence cross-check for the PRL main evidence table.

The submitted PRL argument uses dynesty nested sampling on the official
PTArcade ceffyl density.  This script evaluates the same low-dimensional
integrals with Sobol quasi-Monte-Carlo points in prior-hypercube space and
computes both

  1. direct evidence, log Z = log E_prior[L], and
  2. thermodynamic integration, log Z = int_0^1 <log L>_beta d beta.

This is a narrow cross-check for the 2--3 dimensional main-table models.  It is
not a replacement for a full PTA parallel-tempering chain; it is an independent
evidence estimator for the final official-density likelihood object.
"""

from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import logsumexp
from scipy.stats import qmc


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "code/prl_H4_H5_H6_experiments.py"
OFFICIAL_ROBUSTNESS = ROOT / "results/T2_NG15yr/bayes_factors/prl_official_density_robustness.json"
H4_ROBUSTNESS = ROOT / "results/T2_NG15yr/bayes_factors/prl_H4_env_robustness.json"
OUTDIR = ROOT / "results/T2_NG15yr/bayes_factors"
OUTFILE = OUTDIR / "prl_evidence_ti_qmc_crosscheck.json"
OUT_MD = OUTDIR / "prl_evidence_ti_qmc_crosscheck.md"

M_POWER = 14
N_QMC = 2**M_POWER
SCRAMBLE_SEEDS = [20260421, 20260422, 20260423, 20260424, 20260425]
BETA_GRID = np.concatenate(([0.0], np.geomspace(1.0e-4, 1.0, 80)))

MODEL_ORDER = [
    "baseline_smbhb_ptarcade_bhb_prior",
    "sigw_gauss",
    "sigw_delta",
    "super",
    "smbhb_env_fixed_gamma",
]

MODEL_LABELS = {
    "baseline_smbhb_ptarcade_bhb_prior": "PTArcade BHB-prior baseline",
    "sigw_gauss": "SIGW-Gaussian",
    "sigw_delta": "SIGW-delta",
    "super": "Cosmic superstrings",
    "smbhb_env_fixed_gamma": "Environmental SMBHB, fixed gamma",
}


def load_helper():
    spec = importlib.util.spec_from_file_location("h456", HELPER)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def builder_for(helper, model_key: str):
    if model_key == "baseline_smbhb_ptarcade_bhb_prior":
        return helper.baseline_bhb_prior()
    if model_key == "smbhb_env_fixed_gamma":
        return helper.build_env_model("smbhb_env_fixed_gamma")
    if model_key in {"sigw_gauss", "sigw_delta", "super"}:
        return helper.build_official_spectrum_model(model_key)
    raise ValueError(model_key)


def evaluate_loglikes(helper, model_key: str, scramble_seed: int):
    cpta, parameter_names, prior_transform, model_path = builder_for(helper, model_key)
    sobol = qmc.Sobol(d=cpta.ndim, scramble=True, seed=scramble_seed)
    cube = sobol.random_base2(m=M_POWER)
    loglike = helper.finite_loglike(cpta)
    values = np.empty(N_QMC, dtype=float)
    t0 = time.time()
    for idx, u in enumerate(cube):
        values[idx] = loglike(prior_transform(u))
    runtime = time.time() - t0
    finite = np.isfinite(values)
    if not np.all(finite):
        values = values[finite]
    return {
        "model": model_key,
        "label": MODEL_LABELS[model_key],
        "dim": int(cpta.ndim),
        "parameter_names": parameter_names,
        "ceffyl_param_names": list(cpta.param_names),
        "model_path": str(model_path) if model_path is not None else None,
        "scramble_seed": int(scramble_seed),
        "n_qmc_requested": int(N_QMC),
        "n_finite": int(values.size),
        "finite_fraction": float(values.size / N_QMC),
        "runtime_s": float(runtime),
        "loglike_min": float(values.min()),
        "loglike_max": float(values.max()),
        "loglikes": values,
    }


def evidence_from_loglikes(loglikes: np.ndarray):
    n = loglikes.size
    direct = float(logsumexp(loglikes) - np.log(n))
    beta_means = []
    log_moments = []
    for beta in BETA_GRID:
        weighted = beta * loglikes
        norm = logsumexp(weighted)
        weights = np.exp(weighted - norm)
        beta_means.append(float(np.sum(weights * loglikes)))
        log_moments.append(float(norm - np.log(n)))
    ti = float(trapezoid(beta_means, BETA_GRID))
    return {
        "lnZ_direct_qmc": direct,
        "lnZ_ti_qmc": ti,
        "lnZ_ti_minus_direct": float(ti - direct),
        "beta_grid_size": int(BETA_GRID.size),
        "beta_min_nonzero": float(BETA_GRID[1]),
        "beta_means": beta_means,
        "log_moments": log_moments,
    }


def load_nested_references():
    refs = {}
    if OFFICIAL_ROBUSTNESS.exists():
        data = json.loads(OFFICIAL_ROBUSTNESS.read_text())
        agg = data.get("aggregate", {})
        refs["sigw_gauss"] = agg.get("sigw_gauss_ptarcade", {}).get("lnB_mean")
        refs["sigw_delta"] = agg.get("sigw_delta_ptarcade", {}).get("lnB_mean")
        refs["super"] = agg.get("super_ptarcade", {}).get("lnB_mean")
    if H4_ROBUSTNESS.exists():
        data = json.loads(H4_ROBUSTNESS.read_text())
        refs["smbhb_env_fixed_gamma"] = (
            data.get("aggregate", {}).get("smbhb_env_fixed_gamma", {}).get("lnB_mean")
        )
    return {k: float(v) for k, v in refs.items() if v is not None}


def summarize(seed_results: dict, nested_refs: dict):
    rows = []
    aggregate = {}
    for seed, models in seed_results.items():
        base = models["baseline_smbhb_ptarcade_bhb_prior"]
        for key in MODEL_ORDER:
            r = models[key]
            row = {
                "scramble_seed": int(seed),
                "model": key,
                "label": MODEL_LABELS[key],
                "lnZ_direct_qmc": r["lnZ_direct_qmc"],
                "lnZ_ti_qmc": r["lnZ_ti_qmc"],
                "lnZ_ti_minus_direct": r["lnZ_ti_minus_direct"],
            }
            if key != "baseline_smbhb_ptarcade_bhb_prior":
                row["lnB_direct_qmc_vs_baseline"] = float(
                    r["lnZ_direct_qmc"] - base["lnZ_direct_qmc"]
                )
                row["lnB_ti_qmc_vs_baseline"] = float(
                    r["lnZ_ti_qmc"] - base["lnZ_ti_qmc"]
                )
                if key in nested_refs:
                    row["nested_robust_mean_lnB"] = nested_refs[key]
                    row["direct_residual_vs_nested_mean"] = float(
                        row["lnB_direct_qmc_vs_baseline"] - nested_refs[key]
                    )
                    row["ti_residual_vs_nested_mean"] = float(
                        row["lnB_ti_qmc_vs_baseline"] - nested_refs[key]
                    )
            rows.append(row)

    for key in MODEL_ORDER:
        key_rows = [r for r in rows if r["model"] == key]
        direct_z = np.array([r["lnZ_direct_qmc"] for r in key_rows], dtype=float)
        ti_z = np.array([r["lnZ_ti_qmc"] for r in key_rows], dtype=float)
        item = {
            "label": MODEL_LABELS[key],
            "n_scrambles": len(key_rows),
            "lnZ_direct_mean": float(direct_z.mean()),
            "lnZ_direct_sample_std": float(direct_z.std(ddof=1)),
            "lnZ_ti_mean": float(ti_z.mean()),
            "lnZ_ti_sample_std": float(ti_z.std(ddof=1)),
            "mean_ti_minus_direct": float((ti_z - direct_z).mean()),
            "max_abs_ti_minus_direct": float(np.max(np.abs(ti_z - direct_z))),
        }
        if key != "baseline_smbhb_ptarcade_bhb_prior":
            direct_b = np.array([r["lnB_direct_qmc_vs_baseline"] for r in key_rows], dtype=float)
            ti_b = np.array([r["lnB_ti_qmc_vs_baseline"] for r in key_rows], dtype=float)
            item.update(
                {
                    "lnB_direct_mean": float(direct_b.mean()),
                    "lnB_direct_sample_std": float(direct_b.std(ddof=1)),
                    "lnB_ti_mean": float(ti_b.mean()),
                    "lnB_ti_sample_std": float(ti_b.std(ddof=1)),
                }
            )
            if key in nested_refs:
                item["nested_robust_mean_lnB"] = nested_refs[key]
                item["direct_residual_vs_nested_mean"] = float(direct_b.mean() - nested_refs[key])
                item["ti_residual_vs_nested_mean"] = float(ti_b.mean() - nested_refs[key])
        aggregate[key] = item
    return rows, aggregate


def strip_loglikes(seed_results: dict):
    stripped = {}
    for seed, models in seed_results.items():
        stripped[str(seed)] = {}
        for key, r in models.items():
            stripped[str(seed)][key] = {
                k: v
                for k, v in r.items()
                if k not in {"loglikes", "beta_means", "log_moments"}
            }
    return stripped


def write_md(record: dict):
    lines = [
        "# PRL Evidence Cross-Check: Sobol-QMC Thermodynamic Integration",
        "",
        f"**Generated**: {record['generated']}",
        f"**JSON**: `{OUTFILE.relative_to(ROOT)}`",
        "",
        "## Method",
        "",
        "- Likelihood: official PTArcade `ceffyl` density.",
        "- Models: PTArcade BHB-prior baseline, SIGW-Gaussian, SIGW-delta, cosmic superstrings, and fixed-gamma environmental SMBHB.",
        f"- QMC: `{N_QMC}` Sobol points per scramble, `{len(SCRAMBLE_SEEDS)}` independent scrambles.",
        f"- TI grid: `{len(BETA_GRID)}` beta values from `0` to `1`, with logarithmic spacing above `{BETA_GRID[1]:.1e}`.",
        "- Evidence estimators: direct `log E_prior[L]` and thermodynamic integration `int_0^1 <log L>_beta d beta`.",
        "",
        "## Main Comparison",
        "",
        "| Model | QMC direct lnB | QMC TI lnB | dynesty robust mean | TI residual |",
        "|---|---:|---:|---:|---:|",
    ]
    for key in ["sigw_gauss", "sigw_delta", "super", "smbhb_env_fixed_gamma"]:
        agg = record["aggregate"][key]
        dyn = agg.get("nested_robust_mean_lnB")
        dyn_s = "" if dyn is None else f"`{dyn:+.3f}`"
        resid = agg.get("ti_residual_vs_nested_mean")
        resid_s = "" if resid is None else f"`{resid:+.3f}`"
        lines.append(
            f"| {agg['label']} | `{agg['lnB_direct_mean']:+.3f} +/- {agg['lnB_direct_sample_std']:.3f}` "
            f"| `{agg['lnB_ti_mean']:+.3f} +/- {agg['lnB_ti_sample_std']:.3f}` "
            f"| {dyn_s} | {resid_s} |"
        )
    lines.extend(
        [
            "",
            "## Gate Interpretation",
            "",
            "This cross-check is independent of dynesty and uses a power-posterior",
            "thermodynamic-integration identity evaluated directly in the prior cube.",
            "It is suitable for these 2--3 dimensional official-density integrals.",
            "It is not a full parallel-tempering PTA-chain replacement.",
            "",
            "Diagnostic pass criterion used here: the QMC-TI mean lnB for each main",
            "model must agree with the existing dynesty robust mean to within `0.25` nat,",
            "and the direct-QMC/TI discrepancy must remain below `0.05` nat.",
            "",
            f"**Diagnostic pass**: `{record['diagnostic_pass']}`",
            "",
            "## Remaining Internal-Gate Note",
            "",
            "The project AGENTS rule says Bayes factors must be cross-validated by",
            "thermodynamic integration.  This file supplies a thermodynamic-integration",
            "cross-check for the final official-density likelihood.  If a referee or",
            "internal arbiter specifically requires thermodynamic integration from",
            "parallel-tempered PTA chains, that remains a separate, larger run.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main():
    helper = load_helper()
    nested_refs = load_nested_references()
    seed_results = {}
    print("PRL evidence cross-check: Sobol-QMC + thermodynamic integration")
    print(f"N_QMC={N_QMC}, scrambles={SCRAMBLE_SEEDS}")
    for seed in SCRAMBLE_SEEDS:
        print(f"\n## scramble seed {seed}", flush=True)
        seed_results[seed] = {}
        for key in MODEL_ORDER:
            print(f"running {key}", flush=True)
            raw = evaluate_loglikes(helper, key, seed)
            ev = evidence_from_loglikes(raw.pop("loglikes"))
            raw.update(ev)
            seed_results[seed][key] = raw
            print(
                f"  lnZ direct={raw['lnZ_direct_qmc']:+.3f}; "
                f"lnZ TI={raw['lnZ_ti_qmc']:+.3f}; "
                f"diff={raw['lnZ_ti_minus_direct']:+.3f}; "
                f"runtime={raw['runtime_s']:.1f}s",
                flush=True,
            )

    rows, aggregate = summarize(seed_results, nested_refs)
    max_abs_ti_direct = max(v["max_abs_ti_minus_direct"] for v in aggregate.values())
    residuals = [
        abs(v["ti_residual_vs_nested_mean"])
        for k, v in aggregate.items()
        if k != "baseline_smbhb_ptarcade_bhb_prior" and "ti_residual_vs_nested_mean" in v
    ]
    diagnostic_pass = bool(max_abs_ti_direct <= 0.05 and residuals and max(residuals) <= 0.25)
    record = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "complete": True,
        "task": "PRL main-table independent evidence cross-check",
        "method": {
            "likelihood": "official PTArcade ceffyl density",
            "evidence_estimators": [
                "direct Sobol-QMC prior-space integral",
                "power-posterior thermodynamic integration from same Sobol log-likelihood evaluations",
            ],
            "m_power": M_POWER,
            "n_qmc_per_scramble": N_QMC,
            "scramble_seeds": SCRAMBLE_SEEDS,
            "beta_grid": BETA_GRID.tolist(),
            "models": MODEL_ORDER,
        },
        "nested_references": nested_refs,
        "rows": rows,
        "aggregate": aggregate,
        "diagnostic_pass": diagnostic_pass,
        "diagnostic_pass_criteria": {
            "max_abs_ti_minus_direct_lnZ": 0.05,
            "max_abs_ti_lnb_residual_vs_dynesty_mean": 0.25,
            "note": "This is a low-dimensional official-density TI/QMC cross-check, not a full parallel-tempering PTA-chain TI run.",
        },
        "seed_results": strip_loglikes(seed_results),
    }
    OUTFILE.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n")
    write_md(record)
    print(f"\nsaved: {OUTFILE}")
    print(f"saved: {OUT_MD}")
    print(f"diagnostic_pass={diagnostic_pass}")


if __name__ == "__main__":
    main()
