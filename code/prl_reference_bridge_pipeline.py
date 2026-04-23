#!/usr/bin/env python3
"""P0--P6 bridge pipeline for the PRL/reference-project merge.

This script builds a controlled bridge between:

1. the reference project's public posterior-summary PTA likelihoods; and
2. the current PRL's official NANOGrav ceffyl-density likelihood.

It is intentionally conservative.  Every stage writes machine-readable outputs
and records failed adapters instead of silently dropping them.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import re
import sys
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from ceffyl import Ceffyl
from dynesty import NestedSampler
from dynesty.utils import resample_equal
from enterprise.signals import parameter
from enterprise.signals.parameter import function


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "参考研究课题-PTA引力波"
REF_CODE = REF / "code"
OUT = ROOT / "results/prl_reference_bridge"
OUT.mkdir(parents=True, exist_ok=True)
LOG = OUT / "pipeline.log"

if str(REF_CODE) not in sys.path:
    sys.path.insert(0, str(REF_CODE))

from likelihoods import (  # noqa: E402
    PPTAFreeSpectrumLUT,
    NG15FreeSpectrumLUT,
    EPTAFreeSpectrumLUT,
)
from models import AVAILABLE_MODELS, F_YR  # noqa: E402


DATADIR_OFFICIAL = ROOT / "data/NG15yr/PTArcade_ceffyl_0.2.0/ng15_30f_fs{hd}_ceffyl"
PTARCADE_MODELS = ROOT / "data/NG15yr/PTArcade_models_1.0.0/models_1.0.0"

NSEC_PER_YR = 365.25 * 86400.0
DEFAULT_N_FREQ = 14


PROFILES = {
    "smoke": {
        "nlive": 80,
        "dlogz": 0.5,
        "seed": 20260422,
        "stability": [(80, 20260422), (120, 20260423)],
        "qmc_power": 10,
        "top_n": 4,
    },
    "pilot": {
        "nlive": 250,
        "dlogz": 0.2,
        "seed": 20260422,
        "stability": [(250, 20260422), (250, 20260423), (500, 20260422)],
        "qmc_power": 12,
        "top_n": 5,
    },
    "production": {
        "nlive": 500,
        "dlogz": 0.1,
        "seed": 20260422,
        "stability": [(500, 20260422), (500, 20260423), (1000, 20260422)],
        "qmc_power": 13,
        "top_n": 6,
    },
}


def log(msg: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{stamp}] {msg}"
    print(line, flush=True)
    with LOG.open("a") as f:
        f.write(line + "\n")


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")


def read_json(path: Path):
    return json.loads(path.read_text())


def uniform_prior_from_typename(obj) -> tuple[float, float]:
    typename = getattr(obj, "_typename", "")
    m = re.search(r"Uniform\(pmin=([\-0-9.eE+]+), pmax=([\-0-9.eE+]+)\)", typename)
    if not m:
        raise ValueError(f"only Uniform priors are supported here; got {typename!r}")
    return float(m.group(1)), float(m.group(2))


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


class BaseBridgeModel:
    name = "BaseBridgeModel"
    family = "unknown"
    source_type = "unknown"
    param_names: list[str] = []
    param_defaults: dict[str, float] = {}
    param_priors: dict[str, tuple[float, float]] = {}
    references: list[str] = []

    @property
    def n_params(self) -> int:
        return len(self.param_names)

    def timing_residual_psd(self, f: np.ndarray, **params) -> np.ndarray:
        h2 = self._h_c_squared(f, **params)
        return h2 / (12.0 * np.pi**2 * f**3)

    def omega_gw(self, f: np.ndarray, **params) -> np.ndarray:
        h2 = self._h_c_squared(f, **params)
        # Use the same H0 convention as the reference project for SMBHB controls.
        H0_SI = 2.184e-18
        return (2.0 * np.pi**2 / (3.0 * H0_SI**2)) * f**2 * h2

    def _h_c_squared(self, f: np.ndarray, **params) -> np.ndarray:
        raise NotImplementedError


class CurvedSMBHB(BaseBridgeModel):
    source_type = "astrophysical"
    references = [
        "Sampson et al. 2015, arXiv:1503.02662",
        "NANOGrav 15yr astrophysical interpretations, arXiv:2306.16221",
    ]

    def _base_h2(self, f: np.ndarray, params: dict) -> np.ndarray:
        log10_A = params.get("log10_A", self.param_defaults["log10_A"])
        gamma = params.get("gamma", self.param_defaults.get("gamma", 13.0 / 3.0))
        return 10.0 ** (2.0 * log10_A) * (f / F_YR) ** (3.0 - gamma)


class EnvSMBHB(CurvedSMBHB):
    name = "SMBHB-Env"
    family = "Astro-curved"
    param_names = ["log10_A", "log10_fbend", "kappa"]
    param_defaults = {"log10_A": -14.6, "gamma": 13.0 / 3.0, "log10_fbend": -8.5, "kappa": 2.0}
    param_priors = {"log10_A": (-18.0, -11.0), "log10_fbend": (-9.8, -7.0), "kappa": (0.5, 8.0)}

    def _h_c_squared(self, f: np.ndarray, **params) -> np.ndarray:
        fb = 10.0 ** params.get("log10_fbend", self.param_defaults["log10_fbend"])
        kappa = params.get("kappa", self.param_defaults["kappa"])
        return self._base_h2(f, params) / (1.0 + (fb / f) ** kappa)


class BrokenPLSMBHB(CurvedSMBHB):
    name = "SMBHB-BrokenPL"
    family = "Astro-curved"
    param_names = ["log10_A", "log10_fbend", "delta", "kappa"]
    param_defaults = {
        "log10_A": -14.6,
        "gamma": 13.0 / 3.0,
        "log10_fbend": -8.5,
        "delta": 2.0,
        "kappa": 2.0,
    }
    param_priors = {
        "log10_A": (-18.0, -11.0),
        "log10_fbend": (-10.0, -6.8),
        "delta": (0.0, 6.0),
        "kappa": (0.5, 8.0),
    }

    def _h_c_squared(self, f: np.ndarray, **params) -> np.ndarray:
        fb = 10.0 ** params.get("log10_fbend", self.param_defaults["log10_fbend"])
        delta = params.get("delta", self.param_defaults["delta"])
        kappa = params.get("kappa", self.param_defaults["kappa"])
        return self._base_h2(f, params) * (1.0 + (fb / f) ** kappa) ** (-delta / kappa)


class EccentricSMBHB(CurvedSMBHB):
    name = "SMBHB-Eccentric"
    family = "Astro-curved"
    param_names = ["log10_A", "log10_fe", "beta"]
    param_defaults = {"log10_A": -14.6, "gamma": 13.0 / 3.0, "log10_fe": -8.5, "beta": 1.0}
    param_priors = {"log10_A": (-18.0, -11.0), "log10_fe": (-10.0, -6.8), "beta": (0.0, 3.0)}

    def _h_c_squared(self, f: np.ndarray, **params) -> np.ndarray:
        fe = 10.0 ** params.get("log10_fe", self.param_defaults["log10_fe"])
        beta = params.get("beta", self.param_defaults["beta"])
        return self._base_h2(f, params) * (1.0 + (fe / f) ** 2.0) ** (-beta)


class PTArcadeSpectrumBridge(BaseBridgeModel):
    source_type = "cosmological"

    def __init__(self, name: str, family: str, filename: str):
        self.name = name
        self.family = family
        self.filename = filename
        self.path = PTARCADE_MODELS / filename
        self.module = load_module(self.path, f"bridge_{filename.replace('.', '_')}")
        self.param_names = list(self.module.parameters.keys())
        self.param_priors = {
            key: uniform_prior_from_typename(val)
            for key, val in self.module.parameters.items()
        }
        self.param_defaults = {
            key: 0.5 * (lo + hi) for key, (lo, hi) in self.param_priors.items()
        }
        self.references = [f"PTArcade model file: {self.path}"]

    def timing_residual_psd(self, f: np.ndarray, **params) -> np.ndarray:
        import ptarcade.models_utils as aux

        p = {k: params.get(k, self.param_defaults[k]) for k in self.param_names}
        h2_omega = np.asarray(self.module.spectrum(f, **p), dtype=float)
        h2_omega = np.maximum(h2_omega, 0.0)
        hcf = aux.H_0_Hz / aux.h * np.sqrt(3.0 * h2_omega / 2.0) / (np.pi * f)
        return hcf**2 / (12.0 * np.pi**2 * f**3)

    def omega_gw(self, f: np.ndarray, **params) -> np.ndarray:
        p = {k: params.get(k, self.param_defaults[k]) for k in self.param_names}
        return np.asarray(self.module.spectrum(f, **p), dtype=float)


class ReferenceModelBridge:
    """Thin wrapper around a reference-project GWBSpectrum class."""

    def __init__(self, key: str, family: str):
        self.key = key
        self._model = AVAILABLE_MODELS[key]()
        self.name = key
        self.family = family
        self.source_type = getattr(self._model, "source_type", "unknown")
        self.param_names = list(self._model.param_names)
        self.param_defaults = dict(self._model.param_defaults)
        self.param_priors = dict(self._model.param_priors)
        self.references = list(getattr(self._model, "references", []))

    @property
    def n_params(self) -> int:
        return len(self.param_names)

    def timing_residual_psd(self, f: np.ndarray, **params) -> np.ndarray:
        return self._model.timing_residual_psd(f, **params)

    def omega_gw(self, f: np.ndarray, **params) -> np.ndarray:
        return self._model.omega_gw(f, **params)


def bridge_models() -> dict[str, object]:
    models: dict[str, object] = {
        "SMBHB-PowerLaw": ReferenceModelBridge("SMBHB-PowerLaw", "Astro-simple"),
        "SMBHB-Turnover": ReferenceModelBridge("SMBHB-Turnover", "Astro-simple"),
        "SMBHB-Env": EnvSMBHB(),
        "SMBHB-BrokenPL": BrokenPLSMBHB(),
        "SMBHB-Eccentric": EccentricSMBHB(),
        "PBH-SIGW-Analytic": ReferenceModelBridge("PBH-SIGW-Analytic", "SIGW-like"),
        "SIGW-Gaussian": PTArcadeSpectrumBridge("SIGW-Gaussian", "SIGW-like", "sigw_gauss.py"),
        "SIGW-Delta": PTArcadeSpectrumBridge("SIGW-Delta", "SIGW-like", "sigw_delta.py"),
        "Cosmic-Superstrings": PTArcadeSpectrumBridge("Cosmic-Superstrings", "Cosmo-other", "super.py"),
        "DomainWall": ReferenceModelBridge("DomainWall", "Cosmo-other"),
        "PhaseTransition-Total": ReferenceModelBridge("PhaseTransition-Total", "Cosmo-other"),
    }
    return models


def prior_arrays(model) -> tuple[list[str], np.ndarray, np.ndarray]:
    names = []
    lo = []
    hi = []
    for name in model.param_names:
        default = model.param_defaults.get(name)
        if isinstance(default, str):
            continue
        if name not in model.param_priors:
            raise RuntimeError(f"{model.name}: missing prior for {name}")
        a, b = model.param_priors[name]
        names.append(name)
        lo.append(float(a))
        hi.append(float(b))
    lo_arr = np.asarray(lo, dtype=float)
    hi_arr = np.asarray(hi, dtype=float)
    if np.any(~np.isfinite(lo_arr)) or np.any(~np.isfinite(hi_arr)) or np.any(hi_arr <= lo_arr):
        raise RuntimeError(f"{model.name}: invalid prior ranges")
    return names, lo_arr, hi_arr


def make_prior_transform(model):
    names, lo, hi = prior_arrays(model)
    width = hi - lo

    def pt(u):
        return lo + np.asarray(u) * width

    return names, lo, hi, pt


def params_from_theta(model, names: list[str], theta: Iterable[float]) -> dict[str, float]:
    params = dict(model.param_defaults)
    for name, val in zip(names, theta):
        params[name] = float(val)
    return params


def make_official_loglike(model, n_freq: int = DEFAULT_N_FREQ):
    cpta = Ceffyl.ceffyl(str(DATADIR_OFFICIAL))

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
    order = list(cpta.param_names)

    def ll_from_params(params_dict: dict[str, float]) -> float:
        theta = np.asarray([params_dict[k] for k in order], dtype=float)
        val = cpta.ln_likelihood(theta)
        return float(val) if np.isfinite(val) else -1.0e100

    return ll_from_params, order


@dataclass
class EvidenceResult:
    model: str
    tier: str
    ln_Z: float
    ln_Z_err: float
    n_dim: int
    nlive: int
    dlogz: float
    seed: int
    n_likelihood_calls: int
    wall_seconds: float
    param_names: list[str]
    prior_range: dict[str, list[float]]
    median_params: dict[str, float]
    ci_68: dict[str, list[float]]
    output_dir: str


class BridgeEvidenceProblem:
    def __init__(
        self,
        tier: str,
        model,
        *,
        n_freq: int | dict[str, int] = DEFAULT_N_FREQ,
        ng15_analysis: str = "hd",
    ):
        self.tier = tier
        self.model = model
        self.ng15_analysis = ng15_analysis
        if isinstance(n_freq, dict):
            self.n_freq = n_freq
        else:
            self.n_freq = {"ppta": n_freq, "ng15": n_freq, "epta": n_freq, "official": n_freq}

        self.ppta = None
        self.ng15_local = None
        self.epta = None
        self.official_ll = None
        self.official_order = None

        if tier in ("ng15_local", "local3"):
            self.ng15_local = NG15FreeSpectrumLUT(n_freq=self.n_freq["ng15"], analysis=ng15_analysis)
        if tier in ("ppta", "local3", "hybrid3", "ng15off_ppta"):
            self.ppta = PPTAFreeSpectrumLUT(n_freq=self.n_freq["ppta"])
        if tier in ("epta", "local3", "hybrid3", "ng15off_epta"):
            self.epta = EPTAFreeSpectrumLUT(n_freq=self.n_freq["epta"], dataset="DR2new")
        if tier in ("ng15_official", "hybrid3", "ng15off_ppta", "ng15off_epta"):
            self.official_ll, self.official_order = make_official_loglike(model, self.n_freq["official"])

    def log_likelihood(self, params: dict[str, float]) -> float:
        try:
            if self.tier == "ng15_local":
                return float(self.ng15_local.log_likelihood(self.model, params))
            if self.tier == "ng15_official":
                return float(self.official_ll(params))
            if self.tier == "local3":
                return float(
                    self.ppta.log_likelihood(self.model, params)
                    + self.ng15_local.log_likelihood(self.model, params)
                    + self.epta.log_likelihood(self.model, params)
                )
            if self.tier == "hybrid3":
                return float(
                    self.ppta.log_likelihood(self.model, params)
                    + self.official_ll(params)
                    + self.epta.log_likelihood(self.model, params)
                )
            if self.tier == "ng15off_ppta":
                return float(
                    self.official_ll(params)
                    + self.ppta.log_likelihood(self.model, params)
                )
            if self.tier == "ng15off_epta":
                return float(
                    self.official_ll(params)
                    + self.epta.log_likelihood(self.model, params)
                )
        except Exception:
            return -1.0e100
        raise ValueError(f"unknown tier: {self.tier}")


def run_evidence(
    model,
    tier: str,
    *,
    profile: dict,
    out_dir: Path,
    force: bool = False,
    n_freq: int | dict[str, int] = DEFAULT_N_FREQ,
    seed: int | None = None,
    nlive: int | None = None,
    dlogz: float | None = None,
) -> EvidenceResult:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "ln_Z.json"
    if json_path.exists() and not force:
        d = read_json(json_path)
        return EvidenceResult(**d)

    seed = int(profile["seed"] if seed is None else seed)
    nlive = int(profile["nlive"] if nlive is None else nlive)
    dlogz = float(profile["dlogz"] if dlogz is None else dlogz)
    names, lo, hi, pt = make_prior_transform(model)
    problem = BridgeEvidenceProblem(tier, model, n_freq=n_freq)

    def ll(theta):
        if np.any(~np.isfinite(theta)):
            return -1.0e100
        params = params_from_theta(model, names, theta)
        val = problem.log_likelihood(params)
        return float(val) if np.isfinite(val) else -1.0e100

    log(f"RUN {tier} {model.name} nlive={nlive} dlogz={dlogz} seed={seed}")
    t0 = time.time()
    sampler = NestedSampler(
        ll,
        pt,
        len(names),
        nlive=nlive,
        bound="multi",
        sample="rwalk" if len(names) > 2 else "unif",
        rstate=np.random.default_rng(seed),
    )
    sampler.run_nested(dlogz=dlogz, print_progress=False)
    res = sampler.results
    wall = time.time() - t0
    weights = np.exp(res.logwt - res.logz[-1])
    equal = resample_equal(res.samples, weights)
    med = {names[i]: float(np.median(equal[:, i])) for i in range(len(names))}
    ci68 = {
        names[i]: [float(np.percentile(equal[:, i], 16)), float(np.percentile(equal[:, i], 84))]
        for i in range(len(names))
    }
    np.savez_compressed(out_dir / "posterior_equal_weight.npz", samples=equal, param_names=np.asarray(names))
    result = EvidenceResult(
        model=model.name,
        tier=tier,
        ln_Z=float(res.logz[-1]),
        ln_Z_err=float(res.logzerr[-1]),
        n_dim=len(names),
        nlive=nlive,
        dlogz=dlogz,
        seed=seed,
        n_likelihood_calls=int(np.asarray(res.ncall).sum()),
        wall_seconds=float(wall),
        param_names=names,
        prior_range={names[i]: [float(lo[i]), float(hi[i])] for i in range(len(names))},
        median_params=med,
        ci_68=ci68,
        output_dir=str(out_dir),
    )
    write_json(json_path, result.__dict__)
    log(f"DONE {tier} {model.name}: lnZ={result.ln_Z:+.3f} +/- {result.ln_Z_err:.3f}; wall={wall:.1f}s")
    return result


def csv_write(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})


def stage_p0(args, models: dict[str, object]) -> list[dict]:
    log("P0 bridge manifest")
    rows = []
    f = np.logspace(-10, -6, 64)
    for key, model in models.items():
        row = {
            "model": key,
            "family": getattr(model, "family", ""),
            "source_type": getattr(model, "source_type", ""),
            "n_params": len(model.param_names),
            "param_names": ";".join(model.param_names),
            "prior_ranges": json.dumps(model.param_priors, ensure_ascii=False),
            "local_ppta": "unknown",
            "local_ng15": "unknown",
            "local_epta": "unknown",
            "official_ng15": "unknown",
            "dimension_sanity": "unknown",
            "status": "pending",
            "error": "",
        }
        try:
            psd = model.timing_residual_psd(f, **model.param_defaults)
            omg = model.omega_gw(f, **model.param_defaults)
            finite = np.all(np.isfinite(psd)) and np.all(np.isfinite(omg))
            nonneg = np.all(psd >= 0) and np.all(omg >= 0)
            row["dimension_sanity"] = "pass" if finite and nonneg else "fail"
            # Adapter construction smoke tests.
            PPTAFreeSpectrumLUT(n_freq=2).log_likelihood(model, model.param_defaults)
            NG15FreeSpectrumLUT(n_freq=2, analysis="hd").log_likelihood(model, model.param_defaults)
            EPTAFreeSpectrumLUT(n_freq=2, dataset="DR2new").log_likelihood(model, model.param_defaults)
            official_ll, _ = make_official_loglike(model, n_freq=2)
            official_ll(model.param_defaults)
            row.update({
                "local_ppta": "available",
                "local_ng15": "available_hd",
                "local_epta": "available",
                "official_ng15": "available_hd",
                "status": "pass",
            })
        except Exception as exc:
            row["status"] = "fail"
            row["error"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)

    csv_write(
        OUT / "bridge_model_manifest.csv",
        rows,
        [
            "model",
            "family",
            "source_type",
            "n_params",
            "param_names",
            "prior_ranges",
            "local_ppta",
            "local_ng15",
            "local_epta",
            "official_ng15",
            "dimension_sanity",
            "status",
            "error",
        ],
    )
    passed = [r for r in rows if r["status"] == "pass"]
    md = [
        "# P0 bridge model manifest",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"Passed adapters: {len(passed)} / {len(rows)}",
        "",
        "| Model | Family | n params | status | official NG15 | error |",
        "|---|---|---:|---|---|---|",
    ]
    for r in rows:
        md.append(f"| `{r['model']}` | `{r['family']}` | {r['n_params']} | `{r['status']}` | `{r['official_ng15']}` | {r['error']} |")
    (OUT / "bridge_model_manifest.md").write_text("\n".join(md) + "\n")
    return rows


def load_manifest_passes() -> list[str]:
    path = OUT / "bridge_model_manifest.csv"
    if not path.exists():
        raise FileNotFoundError("P0 manifest is missing")
    with path.open() as f:
        rows = list(csv.DictReader(f))
    return [r["model"] for r in rows if r["status"] == "pass"]


def stage_p1(args, models: dict[str, object], profile: dict) -> None:
    log("P1 NANOGrav calibration anchor")
    keys = load_manifest_passes()
    rows = []
    for key in keys:
        model = models[key]
        try:
            local = run_evidence(
                model,
                "ng15_local",
                profile=profile,
                out_dir=OUT / "P1_ng15_local" / key,
                force=args.force,
            )
            official = run_evidence(
                model,
                "ng15_official",
                profile=profile,
                out_dir=OUT / "P1_ng15_official" / key,
                force=args.force,
            )
            rows.append({
                "model": key,
                "family": getattr(model, "family", ""),
                "lnZ_NG15_local": local.ln_Z,
                "err_NG15_local": local.ln_Z_err,
                "lnZ_NG15_official": official.ln_Z,
                "err_NG15_official": official.ln_Z_err,
                "delta_cal": official.ln_Z - local.ln_Z,
                "delta_cal_err": math.sqrt(local.ln_Z_err**2 + official.ln_Z_err**2),
                "status": "pass",
                "error": "",
            })
        except Exception as exc:
            rows.append({
                "model": key,
                "family": getattr(model, "family", ""),
                "status": "fail",
                "error": f"{type(exc).__name__}: {exc}",
            })
            log(traceback.format_exc())
    csv_write(
        OUT / "ng15_calibration_anchor.csv",
        rows,
        ["model", "family", "lnZ_NG15_local", "err_NG15_local", "lnZ_NG15_official", "err_NG15_official", "delta_cal", "delta_cal_err", "status", "error"],
    )
    write_anchor_md(rows)


def write_anchor_md(rows: list[dict]) -> None:
    ok = [r for r in rows if r.get("status") == "pass"]
    md = [
        "# P1 NANOGrav calibration anchor",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "| Model | Family | lnZ local | lnZ official | delta_cal |",
        "|---|---|---:|---:|---:|",
    ]
    for r in ok:
        md.append(
            f"| `{r['model']}` | `{r['family']}` | {float(r['lnZ_NG15_local']):+.3f} +/- {float(r['err_NG15_local']):.3f} | "
            f"{float(r['lnZ_NG15_official']):+.3f} +/- {float(r['err_NG15_official']):.3f} | "
            f"{float(r['delta_cal']):+.3f} +/- {float(r['delta_cal_err']):.3f} |"
        )
    failed = [r for r in rows if r.get("status") != "pass"]
    if failed:
        md.extend(["", "## Failed rows", ""])
        for r in failed:
            md.append(f"- `{r['model']}`: {r.get('error','')}")
    (OUT / "ng15_calibration_anchor.md").write_text("\n".join(md) + "\n")


def stage_p2(args, models: dict[str, object], profile: dict) -> None:
    log("P2 local3 versus hybrid3 evidence")
    keys = load_manifest_passes()
    for tier in ("local3", "hybrid3"):
        rows = []
        for key in keys:
            model = models[key]
            try:
                res = run_evidence(
                    model,
                    tier,
                    profile=profile,
                    out_dir=OUT / f"P2_{tier}" / key,
                    force=args.force,
                )
                rows.append({
                    "model": key,
                    "family": getattr(model, "family", ""),
                    "ln_Z": res.ln_Z,
                    "ln_Z_err": res.ln_Z_err,
                    "n_dim": res.n_dim,
                    "nlive": res.nlive,
                    "seed": res.seed,
                    "status": "pass",
                    "error": "",
                })
            except Exception as exc:
                rows.append({"model": key, "family": getattr(model, "family", ""), "status": "fail", "error": f"{type(exc).__name__}: {exc}"})
                log(traceback.format_exc())
        rows_sorted = sorted([r for r in rows if r.get("status") == "pass"], key=lambda x: float(x["ln_Z"]), reverse=True)
        best = float(rows_sorted[0]["ln_Z"]) if rows_sorted else np.nan
        for r in rows:
            if r.get("status") == "pass":
                r["rank"] = 1 + rows_sorted.index(r)
                r["delta_to_best"] = float(r["ln_Z"]) - best
        csv_write(OUT / f"{tier}_bridge_lnZ.csv", rows, ["rank", "model", "family", "ln_Z", "ln_Z_err", "delta_to_best", "n_dim", "nlive", "seed", "status", "error"])
    write_ranking_md()


def read_csv_dict(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_ranking_md() -> None:
    md = ["# P2 local3 versus hybrid3 ranking", "", f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}", ""]
    for tier in ("local3", "hybrid3"):
        rows = [r for r in read_csv_dict(OUT / f"{tier}_bridge_lnZ.csv") if r.get("status") == "pass"]
        rows.sort(key=lambda r: int(r["rank"]))
        md.extend([f"## {tier}", "", "| Rank | Model | Family | ln Z | delta to best |", "|---:|---|---|---:|---:|"])
        for r in rows:
            md.append(f"| {r['rank']} | `{r['model']}` | `{r['family']}` | {float(r['ln_Z']):+.3f} +/- {float(r['ln_Z_err']):.3f} | {float(r['delta_to_best']):+.3f} |")
        md.append("")
    (OUT / "local3_vs_hybrid3_ranking.md").write_text("\n".join(md) + "\n")


def stage_p2seq(args, models: dict[str, object], profile: dict) -> None:
    log("P2 sequential anchored bridge ablation")
    keys = load_manifest_passes()
    tier_specs = [
        ("ng15_official", "NG15-off", "ng15_official", OUT / "P1_ng15_official"),
        ("ng15off_ppta", "NG15-off+PPTA-local", "ng15off_ppta", OUT / "P2_sequential" / "ng15off_ppta"),
        ("ng15off_epta", "NG15-off+EPTA-local", "ng15off_epta", OUT / "P2_sequential" / "ng15off_epta"),
        ("hybrid3", "hybrid3", "hybrid3", OUT / "P2_hybrid3"),
    ]
    rows = []
    for tier_key, label, run_tier, base_dir in tier_specs:
        tier_rows = []
        for key in keys:
            model = models[key]
            try:
                res = run_evidence(
                    model,
                    run_tier,
                    profile=profile,
                    out_dir=base_dir / key,
                    force=args.force and tier_key not in ("ng15_official", "hybrid3"),
                )
                row = {
                    "tier": tier_key,
                    "label": label,
                    "model": key,
                    "family": getattr(model, "family", ""),
                    "ln_Z": res.ln_Z,
                    "ln_Z_err": res.ln_Z_err,
                    "n_dim": res.n_dim,
                    "nlive": res.nlive,
                    "seed": res.seed,
                    "status": "pass",
                    "error": "",
                }
            except Exception as exc:
                row = {
                    "tier": tier_key,
                    "label": label,
                    "model": key,
                    "family": getattr(model, "family", ""),
                    "status": "fail",
                    "error": f"{type(exc).__name__}: {exc}",
                }
                log(traceback.format_exc())
            tier_rows.append(row)
        rows.extend(tier_rows)

    for tier_key, _, _, _ in tier_specs:
        good = [r for r in rows if r["tier"] == tier_key and r.get("status") == "pass"]
        good.sort(key=lambda r: float(r["ln_Z"]), reverse=True)
        best = float(good[0]["ln_Z"]) if good else np.nan
        rank = {r["model"]: i + 1 for i, r in enumerate(good)}
        for r in rows:
            if r["tier"] == tier_key and r.get("status") == "pass":
                r["rank"] = rank[r["model"]]
                r["delta_to_best"] = float(r["ln_Z"]) - best

    csv_write(
        OUT / "sequential_bridge_ablation.csv",
        rows,
        ["tier", "label", "rank", "model", "family", "ln_Z", "ln_Z_err", "delta_to_best", "n_dim", "nlive", "seed", "status", "error"],
    )
    write_sequential_bridge_md(tier_specs)


def write_sequential_bridge_md(tier_specs: list[tuple[str, str, str, Path]]) -> None:
    rows = read_csv_dict(OUT / "sequential_bridge_ablation.csv")
    md = [
        "# P2b sequential anchored bridge ablation",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "The sequence is an anchored source-ranking stress test, not a fully calibrated multi-PTA Bayes-factor scale.",
        "",
    ]
    for tier_key, label, _, _ in tier_specs:
        good = [r for r in rows if r["tier"] == tier_key and r.get("status") == "pass"]
        good.sort(key=lambda r: int(r["rank"]))
        md.extend([f"## {label}", "", "| Rank | Model | Family | ln Z | delta to best |", "|---:|---|---|---:|---:|"])
        for r in good:
            md.append(f"| {r['rank']} | `{r['model']}` | `{r['family']}` | {float(r['ln_Z']):+.3f} +/- {float(r['ln_Z_err']):.3f} | {float(r['delta_to_best']):+.3f} |")
        md.append("")
    leaders = []
    for tier_key, label, _, _ in tier_specs:
        good = [r for r in rows if r["tier"] == tier_key and r.get("status") == "pass"]
        if not good:
            continue
        best = min(good, key=lambda r: int(r["rank"]))
        leaders.append(f"- {label}: `{best['model']}` ({float(best['ln_Z']):+.3f})")
    md.extend(["## Leaders", "", *leaders, ""])
    (OUT / "sequential_bridge_ablation.md").write_text("\n".join(md) + "\n")


def logsumexp_family(rows: list[dict], family: str, leave_out: str | None = None) -> tuple[float, int]:
    vals = []
    for r in rows:
        if r["family"] != family or r["model"] == leave_out:
            continue
        vals.append(float(r["ln_Z"]))
    if not vals:
        return np.nan, 0
    # Equal family prior mass; equal within-family weights.
    return float(logsumexp(vals) - np.log(len(vals))), len(vals)


def stage_p3(args) -> None:
    log("P3 family evidence")
    rows_out = []
    families = ["Astro-simple", "Astro-curved", "SIGW-like", "Cosmo-other"]
    for tier in ("local3", "hybrid3"):
        rows = [r for r in read_csv_dict(OUT / f"{tier}_bridge_lnZ.csv") if r.get("status") == "pass"]
        for fam in families:
            z, n = logsumexp_family(rows, fam)
            rows_out.append({"tier": tier, "family": fam, "mode": "all", "leave_out": "", "n_models": n, "ln_Z_family": z})
            fam_rows = [r for r in rows if r["family"] == fam]
            if fam_rows:
                leader = max(fam_rows, key=lambda r: float(r["ln_Z"]))["model"]
                z2, n2 = logsumexp_family(rows, fam, leave_out=leader)
                rows_out.append({"tier": tier, "family": fam, "mode": "leave_one_leader_out", "leave_out": leader, "n_models": n2, "ln_Z_family": z2})
    csv_write(OUT / "family_evidence.csv", rows_out, ["tier", "family", "mode", "leave_out", "n_models", "ln_Z_family"])
    md = ["# P3 family evidence", "", f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}", "", "| Tier | Family | Mode | leave out | n | ln Z_F |", "|---|---|---|---|---:|---:|"]
    for r in rows_out:
        val = "" if not np.isfinite(float(r["ln_Z_family"])) else f"{float(r['ln_Z_family']):+.3f}"
        md.append(f"| `{r['tier']}` | `{r['family']}` | `{r['mode']}` | `{r['leave_out']}` | {r['n_models']} | {val} |")
    (OUT / "family_evidence_sensitivity.md").write_text("\n".join(md) + "\n")


def fit_effective_gamma(model, params: dict, f_bins: np.ndarray, weights: np.ndarray | None = None) -> tuple[float, float, float]:
    log_rho = 0.5 * np.log10(np.maximum(model.timing_residual_psd(f_bins, **params), 1e-300))
    if weights is None:
        weights = np.ones_like(f_bins)
    logf = np.log10(f_bins / F_YR)

    def loss(x):
        log10_A, gamma = x
        pred_h2 = 2.0 * log10_A + (3.0 - gamma) * logf
        pred_log_rho = 0.5 * np.log10(np.maximum((10.0 ** pred_h2) / (12.0 * np.pi**2 * f_bins**3), 1e-300))
        return float(np.sum(weights * (log_rho - pred_log_rho) ** 2))

    res = minimize(loss, x0=np.array([-14.5, 13.0 / 3.0]), bounds=[(-18, -11), (0, 7)], method="Nelder-Mead")
    return float(res.x[0]), float(res.x[1]), float(res.fun)


def stage_p4(args, models: dict[str, object]) -> None:
    log("P4 gamma projection")
    rows = [r for r in read_csv_dict(OUT / "hybrid3_bridge_lnZ.csv") if r.get("status") == "pass"]
    curved = [r["model"] for r in rows if r["family"] == "Astro-curved"]
    # Use top curved models first.
    curved = [r["model"] for r in sorted([r for r in rows if r["family"] == "Astro-curved"], key=lambda x: float(x["ln_Z"]), reverse=True)]
    like_bins = {
        "PPTA-DR3": PPTAFreeSpectrumLUT(n_freq=DEFAULT_N_FREQ).f_bins,
        "NANOGrav-15yr": NG15FreeSpectrumLUT(n_freq=DEFAULT_N_FREQ, analysis="hd").f_bins,
        "EPTA-DR2new": EPTAFreeSpectrumLUT(n_freq=DEFAULT_N_FREQ, dataset="DR2new").f_bins,
    }
    out_rows = []
    for key in curved:
        model = models[key]
        sample_path = OUT / "P2_hybrid3" / key / "posterior_equal_weight.npz"
        if not sample_path.exists():
            continue
        data = np.load(sample_path, allow_pickle=True)
        samples = data["samples"]
        pnames = [str(x) for x in data["param_names"]]
        idx = np.linspace(0, samples.shape[0] - 1, min(512, samples.shape[0]), dtype=int)
        for dataset, f_bins in like_bins.items():
            gammas = []
            amps = []
            losses = []
            for ii in idx:
                params = dict(model.param_defaults)
                for j, pn in enumerate(pnames):
                    params[pn] = float(samples[ii, j])
                amp, gamma, loss = fit_effective_gamma(model, params, f_bins)
                gammas.append(gamma)
                amps.append(amp)
                losses.append(loss)
            out_rows.append({
                "model": key,
                "dataset": dataset,
                "n_projected_samples": len(gammas),
                "gamma_eff_median": float(np.median(gammas)),
                "gamma_eff_p16": float(np.percentile(gammas, 16)),
                "gamma_eff_p84": float(np.percentile(gammas, 84)),
                "log10_A_eff_median": float(np.median(amps)),
                "projection_loss_median": float(np.median(losses)),
            })
    csv_write(OUT / "gamma_projection.csv", out_rows, ["model", "dataset", "n_projected_samples", "gamma_eff_median", "gamma_eff_p16", "gamma_eff_p84", "log10_A_eff_median", "projection_loss_median"])
    md = ["# P4 gamma projection", "", f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}", "", "| Model | Dataset | gamma_eff median | 68% interval |", "|---|---|---:|---:|"]
    for r in out_rows:
        md.append(f"| `{r['model']}` | `{r['dataset']}` | {float(r['gamma_eff_median']):.3f} | [{float(r['gamma_eff_p16']):.3f}, {float(r['gamma_eff_p84']):.3f}] |")
    (OUT / "gamma_projection.md").write_text("\n".join(md) + "\n")


def nfreq_for_cut(freqs: np.ndarray, fmax_nhz: float) -> int:
    return max(1, int(np.sum(freqs <= fmax_nhz * 1e-9)))


def stage_p5(args, models: dict[str, object], profile: dict) -> None:
    log("P5 low-frequency discrimination scan")
    ppta_ref = PPTAFreeSpectrumLUT(n_freq=DEFAULT_N_FREQ)
    ng15_ref = NG15FreeSpectrumLUT(n_freq=DEFAULT_N_FREQ, analysis="hd")
    epta_ref = EPTAFreeSpectrumLUT(n_freq=DEFAULT_N_FREQ, dataset="DR2new")
    p2_rows = [r for r in read_csv_dict(OUT / "hybrid3_bridge_lnZ.csv") if r.get("status") == "pass"]
    p2_rows.sort(key=lambda r: float(r["ln_Z"]), reverse=True)
    top = [r["model"] for r in p2_rows[: int(profile["top_n"])]]
    cuts = [5.0, 8.0, 12.0, 20.0, 1e9]
    out_rows = []
    for cut in cuts:
        nfd = {
            "ppta": nfreq_for_cut(ppta_ref.f_bins, cut),
            "ng15": nfreq_for_cut(ng15_ref.f_bins, cut),
            "epta": nfreq_for_cut(epta_ref.f_bins, cut),
            "official": nfreq_for_cut(ng15_ref.f_bins, cut),
        }
        label = "all" if cut > 1e8 else f"le_{cut:g}_nHz"
        for key in top:
            try:
                res = run_evidence(
                    models[key],
                    "hybrid3",
                    profile=profile,
                    out_dir=OUT / "P5_frequency_cut" / label / key,
                    force=args.force,
                    n_freq=nfd,
                    nlive=max(80, int(profile["nlive"] // 2)),
                    dlogz=max(float(profile["dlogz"]), 0.2),
                )
                out_rows.append({
                    "fmax_nHz": label,
                    "model": key,
                    "family": getattr(models[key], "family", ""),
                    "n_ppta": nfd["ppta"],
                    "n_ng15": nfd["ng15"],
                    "n_epta": nfd["epta"],
                    "ln_Z": res.ln_Z,
                    "ln_Z_err": res.ln_Z_err,
                    "status": "pass",
                    "error": "",
                })
            except Exception as exc:
                out_rows.append({"fmax_nHz": label, "model": key, "family": getattr(models[key], "family", ""), "status": "fail", "error": f"{type(exc).__name__}: {exc}"})
                log(traceback.format_exc())
    # Add delta-to-best within each cut.
    for label in sorted({r["fmax_nHz"] for r in out_rows}):
        good = [r for r in out_rows if r["fmax_nHz"] == label and r.get("status") == "pass"]
        if not good:
            continue
        best = max(float(r["ln_Z"]) for r in good)
        for r in good:
            r["delta_to_best"] = float(r["ln_Z"]) - best
    csv_write(OUT / "frequency_cut_evidence.csv", out_rows, ["fmax_nHz", "model", "family", "n_ppta", "n_ng15", "n_epta", "ln_Z", "ln_Z_err", "delta_to_best", "status", "error"])
    md = ["# P5 frequency-cut evidence", "", f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}", "", "| fmax | Model | ln Z | delta to best |", "|---|---|---:|---:|"]
    for r in out_rows:
        if r.get("status") == "pass":
            md.append(f"| `{r['fmax_nHz']}` | `{r['model']}` | {float(r['ln_Z']):+.3f} +/- {float(r['ln_Z_err']):.3f} | {float(r.get('delta_to_best', 0)):+.3f} |")
    (OUT / "frequency_cut_evidence.md").write_text("\n".join(md) + "\n")


def qmc_ti(problem: BridgeEvidenceProblem, model, power: int, seed: int) -> dict:
    names, _, _, pt = make_prior_transform(model)
    try:
        from scipy.stats import qmc
    except Exception as exc:
        raise RuntimeError(f"QMC unavailable: {exc}")
    sampler = qmc.Sobol(d=len(names), scramble=True, seed=seed)
    cube = sampler.random_base2(m=power)
    lls = np.empty(cube.shape[0], dtype=float)
    for i, u in enumerate(cube):
        theta = pt(u)
        params = params_from_theta(model, names, theta)
        lls[i] = problem.log_likelihood(params)
    finite = lls[np.isfinite(lls)]
    direct = float(logsumexp(finite) - np.log(finite.size))
    betas = np.concatenate(([0.0], np.geomspace(1e-4, 1.0, 64)))
    means = []
    for beta in betas:
        lw = beta * finite
        norm = logsumexp(lw)
        w = np.exp(lw - norm)
        means.append(float(np.sum(w * finite)))
    ti = float(np.trapz(means, betas))
    return {"n_qmc": int(finite.size), "lnZ_qmc_direct": direct, "lnZ_qmc_ti": ti, "loglike_min": float(finite.min()), "loglike_max": float(finite.max())}


def stage_p6(args, models: dict[str, object], profile: dict) -> None:
    log("P6 robustness suite")
    rows = [r for r in read_csv_dict(OUT / "hybrid3_bridge_lnZ.csv") if r.get("status") == "pass"]
    rows.sort(key=lambda r: float(r["ln_Z"]), reverse=True)
    top = [r["model"] for r in rows[: min(4, len(rows))]]
    robust_rows = []
    for key in top:
        model = models[key]
        for nlive, seed in profile["stability"]:
            try:
                res = run_evidence(
                    model,
                    "hybrid3",
                    profile=profile,
                    out_dir=OUT / "P6_stability" / f"nlive{nlive}_seed{seed}" / key,
                    force=args.force,
                    nlive=nlive,
                    seed=seed,
                )
                robust_rows.append({"check": "seed_nlive", "model": key, "nlive": nlive, "seed": seed, "ln_Z": res.ln_Z, "ln_Z_err": res.ln_Z_err, "status": "pass", "error": ""})
            except Exception as exc:
                robust_rows.append({"check": "seed_nlive", "model": key, "nlive": nlive, "seed": seed, "status": "fail", "error": f"{type(exc).__name__}: {exc}"})
    csv_write(OUT / "robustness_budget.csv", robust_rows, ["check", "model", "nlive", "seed", "ln_Z", "ln_Z_err", "status", "error"])

    ti_rows = []
    for key in top[:3]:
        model = models[key]
        try:
            problem = BridgeEvidenceProblem("hybrid3", model)
            out = qmc_ti(problem, model, int(profile["qmc_power"]), seed=20260422)
            out.update({"model": key, "tier": "hybrid3", "status": "pass", "error": ""})
            ti_rows.append(out)
        except Exception as exc:
            ti_rows.append({"model": key, "tier": "hybrid3", "status": "fail", "error": f"{type(exc).__name__}: {exc}"})
            log(traceback.format_exc())
    csv_write(OUT / "ti_crosscheck_top_models.csv", ti_rows, ["model", "tier", "n_qmc", "lnZ_qmc_direct", "lnZ_qmc_ti", "loglike_min", "loglike_max", "status", "error"])
    write_run_summary()


def write_run_summary() -> None:
    lines = [
        "# P0--P6 run summary",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Outputs",
        "",
    ]
    for name in [
        "bridge_model_manifest.csv",
        "ng15_calibration_anchor.csv",
        "local3_bridge_lnZ.csv",
        "hybrid3_bridge_lnZ.csv",
        "family_evidence.csv",
        "gamma_projection.csv",
        "frequency_cut_evidence.csv",
        "robustness_budget.csv",
        "ti_crosscheck_top_models.csv",
        "sequential_bridge_ablation.csv",
    ]:
        p = OUT / name
        lines.append(f"- `{p.relative_to(ROOT)}`: {'exists' if p.exists() else 'missing'}")
    lines.extend(["", "## Next discussion gate", "", "Review whether curved-SMBHB controls enter the hybrid3 top tier and whether family-level evidence remains order-unity separated."])
    (OUT / "P0_P6_run_summary.md").write_text("\n".join(lines) + "\n")


def run_selected(stage: str, args) -> None:
    profile = PROFILES[args.profile]
    models = bridge_models()
    if stage in ("p0", "all"):
        stage_p0(args, models)
    if stage in ("p1", "all"):
        stage_p1(args, models, profile)
    if stage in ("p2", "all"):
        stage_p2(args, models, profile)
    if stage in ("p2seq", "all"):
        stage_p2seq(args, models, profile)
    if stage in ("p3", "all"):
        stage_p3(args)
    if stage in ("p4", "all"):
        stage_p4(args, models)
    if stage in ("p5", "all"):
        stage_p5(args, models, profile)
    if stage in ("p6", "all"):
        stage_p6(args, models, profile)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["p0", "p1", "p2", "p2seq", "p3", "p4", "p5", "p6", "all"])
    parser.add_argument("--profile", choices=sorted(PROFILES), default="production")
    parser.add_argument("--force", action="store_true", help="rerun even if output files exist")
    args = parser.parse_args()
    log(f"START stage={args.stage} profile={args.profile} force={args.force}")
    run_selected(args.stage, args)
    log(f"FINISH stage={args.stage} profile={args.profile}")


if __name__ == "__main__":
    main()
