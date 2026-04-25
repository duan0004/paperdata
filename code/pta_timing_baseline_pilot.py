#!/usr/bin/env python3
"""P2 timing-level baseline pilot for public PTA products.

This is a mechanical ENTERPRISE gate, not a Bayes-factor calculation.  It
loads a small public-data subset, builds a timing-marginalized common
power-law model with fixed public white-noise parameters where available, and
checks that the PTA object has finite prior and likelihood values.

The script does not modify timing-model files and does not implement timing
residuals by hand.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = ROOT / "code"
OUT_DIR = ROOT / "results" / "5pta_timing"

if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import pta_loader_smoke as loaders  # noqa: E402


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def update_noise_dict(target: dict, path: Path) -> None:
    if path.exists():
        target.update(read_json(path))


def load_noise_dicts(arrays: list[str], pulsar_names: list[str], epta_variant: str) -> tuple[dict, dict]:
    noise: dict = {}
    provenance: dict[str, list[str] | str] = {}

    if "NG15" in arrays or "NG15_TIMING" in arrays:
        ng15_noise = ROOT / "data" / "NG15yr" / "tutorials" / "data" / "15yr_wn_dict.json"
        update_noise_dict(noise, ng15_noise)
        provenance["NG15"] = rel(ng15_noise)

    if "EPTA" in arrays:
        epta_root = ROOT / "data" / "EPTA_DR2" / "epta-dr2" / "EPTA-DR2" / "noisefiles" / epta_variant
        used = []
        for psr in pulsar_names:
            path = epta_root / f"{psr}_noise.json"
            if path.exists():
                update_noise_dict(noise, path)
                used.append(rel(path))
        provenance["EPTA"] = used

    if "PPTA" in arrays:
        ppta_root = ROOT / "data" / "PPTA_DR3" / "ppta_dr3" / "toas_and_parameters" / "noisefiles"
        used = []
        for psr in pulsar_names:
            matches = sorted(ppta_root.glob(f"{psr}*_noise.json"))
            for path in matches:
                update_noise_dict(noise, path)
                used.append(rel(path))
        provenance["PPTA"] = used

    return noise, provenance


def _take_unique(
    rows: list,
    max_per_array: int,
    used_names: set[str],
    skipped: list[dict],
    array: str,
    excluded_names: set[str],
):
    selected = []
    for row in rows:
        name = row[0]
        if name in excluded_names:
            skipped.append({"array": array, "name": name, "reason": "excluded by --exclude-pulsars"})
            continue
        if name in used_names:
            skipped.append({"array": array, "name": name, "reason": "duplicate pulsar name in pilot namespace"})
            continue
        selected.append(row)
        used_names.add(name)
        if max_per_array > 0 and len(selected) >= max_per_array:
            break
    return selected


def load_pulsars(
    arrays: list[str],
    max_per_array: int,
    epta_variant: str,
    ephem: str,
    excluded_names: set[str],
) -> tuple[list, list[dict], list[dict]]:
    from enterprise.pulsar import Pulsar
    from enterprise_extensions.load_feathers import load_feathers_from_folder

    psrs = []
    meta = []
    skipped = []
    used_names: set[str] = set()
    ephem_arg = None if ephem == "FROM_PAR" else ephem

    if "NG15" in arrays:
        feather_dir = ROOT / "data" / "NG15yr" / "tutorials" / "data" / "feathers"
        ng15 = load_feathers_from_folder(str(feather_dir))
        rows = [(p.name, p) for p in ng15]
        for _, p in _take_unique(rows, max_per_array, used_names, skipped, "NG15", excluded_names):
            psrs.append(p)
            meta.append({"array": "NG15", "name": p.name, "source": "feather", "parfile": "", "timfile": ""})

    if "NG15_TIMING" in arrays:
        for psr, par, tim in _take_unique(loaders.ng15_timing_pairs(), max_per_array, used_names, skipped, "NG15_TIMING", excluded_names):
            pobj = Pulsar(str(par), str(tim), ephem=ephem_arg)
            psrs.append(pobj)
            meta.append({"array": "NG15_TIMING", "name": pobj.name, "source": "par_tim", "parfile": rel(par), "timfile": rel(tim)})

    if "EPTA" in arrays:
        for psr, par, tim in _take_unique(loaders.epta_pairs(epta_variant), max_per_array, used_names, skipped, "EPTA", excluded_names):
            pobj = Pulsar(str(par), str(tim), ephem=ephem_arg)
            psrs.append(pobj)
            meta.append({"array": "EPTA", "name": pobj.name, "source": "par_tim", "parfile": rel(par), "timfile": rel(tim)})

    if "PPTA" in arrays:
        for psr, par, tim in _take_unique(loaders.ppta_pairs(), max_per_array, used_names, skipped, "PPTA", excluded_names):
            pobj = Pulsar(str(par), str(tim), ephem=ephem_arg)
            psrs.append(pobj)
            meta.append({"array": "PPTA", "name": pobj.name, "source": "par_tim", "parfile": rel(par), "timfile": rel(tim)})

    return psrs, meta, skipped


def ppta_basis_ecorr_selection(flags):
    """PPTA DR3 ECORR keys are grouped by coarse UWL band."""
    vals = np.asarray(flags.get("B", []), dtype=str)
    out = {}
    for band in np.unique(vals):
        if not band or band == "None":
            continue
        if band.startswith("uwl_"):
            key = "basis_ecorr_" + band.replace("uwl_", "") + "_uwl"
        else:
            key = "basis_ecorr_" + band
        out[key] = vals == band
    return out


def canonical_ppta_noise_group(psr_name: str, group: str, noise: dict) -> str:
    """Map PPTA TIM groups to staged noise-product groups without inventing values."""
    if (
        f"{psr_name}_{group}_efac" in noise
        or f"{psr_name}_{group}_log10_tnequad" in noise
        or f"{psr_name}_{group}_log10_t2equad" in noise
    ):
        return group

    pdfb_merged = re.sub(r"^PDFB\d+_", "PDFB_", group)
    if pdfb_merged != group and (
        f"{psr_name}_{pdfb_merged}_efac" in noise
        or f"{psr_name}_{pdfb_merged}_log10_tnequad" in noise
        or f"{psr_name}_{pdfb_merged}_log10_t2equad" in noise
    ):
        return pdfb_merged

    return group


def ppta_white_noise_selection(psr_name: str, noise: dict):
    """Use PPTA `-group` flags, with documented aliases to staged noise keys."""

    def selection(backend_flags, flags):
        groups = np.asarray(flags.get("group", backend_flags), dtype=str)
        out: dict[str, np.ndarray] = {}
        for group in np.unique(groups):
            if not group or group == "None":
                continue
            key = canonical_ppta_noise_group(psr_name, group, noise)
            mask = groups == group
            if key in out:
                out[key] = out[key] | mask
            else:
                out[key] = mask
        return out

    return selection


def ppta_noise_aliases(psrs: list, meta: list[dict], noise: dict) -> list[dict]:
    aliases = []
    for p, item in zip(psrs, meta):
        if item["array"] != "PPTA":
            continue
        groups = np.asarray(p.flags.get("group", p.backend_flags), dtype=str)
        for group in sorted(set(groups)):
            if not group or group == "None":
                continue
            mapped = canonical_ppta_noise_group(p.name, group, noise)
            if mapped != group:
                aliases.append({"pulsar": p.name, "from_group": group, "to_noise_group": mapped, "n_toas": int(np.sum(groups == group))})
    return aliases


def build_pta(psrs: list, meta: list[dict], components_gw: int, components_rn: int, include_red_noise: bool, noise: dict):
    from enterprise.signals import gp_signals, parameter, selections, signal_base, utils, white_signals

    tmin = min(float(np.min(p.toas)) for p in psrs)
    tmax = max(float(np.max(p.toas)) for p in psrs)
    tspan = tmax - tmin

    by_backend = selections.Selection(selections.by_backend)
    ppta_ecorr = selections.Selection(ppta_basis_ecorr_selection)

    log10_A_gw = parameter.Uniform(-18, -13)("log10_A_gw")
    gamma_gw = parameter.Uniform(0, 7)("gamma_gw")
    common_psd = utils.powerlaw(log10_A=log10_A_gw, gamma=gamma_gw)
    gw = gp_signals.FourierBasisCommonGP(
        spectrum=common_psd,
        orf=utils.hd_orf(),
        components=components_gw,
        Tspan=tspan,
        name="gw",
    )

    tm = gp_signals.MarginalizingTimingModel(use_svd=True)

    signals = []
    for item in meta:
        array = item["array"]
        if array in {"NG15", "NG15_TIMING"}:
            mn = white_signals.MeasurementNoise(
                efac=parameter.Constant(),
                log10_t2equad=parameter.Constant(),
                selection=by_backend,
            )
            ec = white_signals.EcorrKernelNoise(log10_ecorr=parameter.Constant(), selection=by_backend)
            sig = tm + mn + ec
        elif array == "EPTA":
            mn = white_signals.MeasurementNoise(efac=parameter.Constant(), selection=by_backend)
            eq = white_signals.TNEquadNoise(log10_tnequad=parameter.Constant(), selection=by_backend)
            sig = tm + mn + eq
        elif array == "PPTA":
            ppta_group = selections.Selection(ppta_white_noise_selection(item["name"], noise))
            mn = white_signals.MeasurementNoise(efac=parameter.Constant(), selection=ppta_group)
            eq = white_signals.TNEquadNoise(log10_tnequad=parameter.Constant(), selection=ppta_group)
            ec = white_signals.EcorrKernelNoise(log10_ecorr=parameter.Constant(), selection=ppta_ecorr)
            sig = tm + mn + eq + ec
        else:
            raise ValueError(f"Unsupported array: {array}")

        if include_red_noise:
            log10_A_rn = parameter.Uniform(-20, -11)
            gamma_rn = parameter.Uniform(0, 7)
            red_psd = utils.powerlaw(log10_A=log10_A_rn, gamma=gamma_rn)
            rn = gp_signals.FourierBasisGP(spectrum=red_psd, components=components_rn, Tspan=tspan)
            sig = sig + rn

        signals.append(sig + gw)

    pta = signal_base.PTA([sig(p) for sig, p in zip(signals, psrs)])
    return pta, tspan


def make_initial_point(pta, defaults: dict, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    x0 = np.asarray([p.sample() for p in pta.params], dtype=float).flatten()
    for i, name in enumerate(pta.param_names):
        if name == "log10_A_gw":
            x0[i] = -14.2
        elif name == "gamma_gw":
            x0[i] = 3.25
        elif name.endswith("_red_noise_log10_A") or name.endswith("_rn_log10_A"):
            candidate = defaults.get(name)
            if isinstance(candidate, (int, float)) and -20 <= float(candidate) <= -11:
                x0[i] = float(candidate)
            else:
                x0[i] = -15.0 + 0.2 * rng.normal()
        elif name.endswith("_red_noise_gamma") or name.endswith("_rn_gamma"):
            candidate = defaults.get(name)
            if isinstance(candidate, (int, float)) and 0 <= float(candidate) <= 7:
                x0[i] = float(candidate)
            else:
                x0[i] = 4.0 + 0.2 * rng.normal()
        elif name in defaults and isinstance(defaults[name], (int, float)) and math.isfinite(float(defaults[name])):
            x0[i] = float(defaults[name])
    return x0


def white_noise_params(param_names: list[str]) -> list[str]:
    suffixes = (
        "_efac",
        "_log10_tnequad",
        "_log10_t2equad",
        "_log10_ecorr",
    )
    return [name for name in param_names if name.endswith(suffixes) or "_basis_ecorr_" in name]


def write_markdown(path: Path, result: dict) -> None:
    lines = [
        "# P2 Timing-Level Baseline Pilot",
        "",
        f"Generated: {result['generated']}",
        "",
        "This is a mechanical ENTERPRISE likelihood gate, not a Bayes-factor or source-identification result.",
        "",
        "## Status",
        "",
        f"- Status: `{result['status']}`",
        f"- Output tag: `{result['output_tag']}`",
        f"- Arrays: `{', '.join(result['arrays'])}`",
        f"- Pulsars loaded: `{result['n_pulsars']}`",
        f"- Duplicate pulsars skipped in pilot namespace: `{len(result.get('skipped_pulsars', []))}`",
        f"- GWB Fourier components: `{result['components_gw']}`",
        f"- Red noise included: `{result['include_red_noise']}`",
        f"- Excluded pulsars: `{', '.join(result.get('exclude_pulsars', [])) or 'none'}`",
        f"- Tspan years: `{result.get('tspan_years', 'NA')}`",
        f"- PTA ndim: `{result.get('ndim', 'NA')}`",
        f"- ln prior at pilot point: `{result.get('lnprior')}`",
        f"- ln likelihood at pilot point: `{result.get('lnlikelihood')}`",
        f"- Missing fixed white-noise parameters: `{len(result.get('missing_white_noise_parameters', []))}`",
        f"- PPTA group aliases applied: `{len(result.get('ppta_noise_aliases', []))}`",
        "",
        "## Public Inputs",
        "",
        "- NANOGrav 15-year data products must be cited with arXiv:2306.16213 and arXiv:2306.16214 when used.",
        f"- Noise parameters loaded: `{result['n_noise_parameters']}` entries.",
        "",
        "## Pulsars",
        "",
        "| Array | Pulsar | source |",
        "|---|---|---|",
    ]
    for p in result["pulsars"]:
        lines.append(f"| {p['array']} | `{p['name']}` | `{p['source']}` |")
    if result.get("missing_white_noise_parameters"):
        lines.extend(["", "## Missing Fixed White-Noise Parameters", ""])
        for name in result["missing_white_noise_parameters"]:
            lines.append(f"- `{name}`")
    if result.get("ppta_noise_aliases"):
        lines.extend(["", "## PPTA Noise Group Aliases", "", "| Pulsar | TIM group | noise group | TOAs |", "|---|---|---|---:|"])
        for row in result["ppta_noise_aliases"]:
            lines.append(f"| `{row['pulsar']}` | `{row['from_group']}` | `{row['to_noise_group']}` | {row['n_toas']} |")
    if result.get("error"):
        lines.extend(["", "## Error", "", "```text", result["error"].rstrip(), "```"])
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arrays", nargs="+", default=["NG15", "EPTA", "PPTA"], choices=["NG15", "NG15_TIMING", "EPTA", "PPTA"])
    parser.add_argument("--max-per-array", type=int, default=2)
    parser.add_argument("--components-gw", type=int, default=4)
    parser.add_argument("--components-rn", type=int, default=5)
    parser.add_argument("--include-red-noise", action="store_true")
    parser.add_argument("--epta-variant", default="DR2new+", choices=["DR2new+", "DR2new", "DR2full+", "DR2full"])
    parser.add_argument("--ephem", default="FROM_PAR", help="Use FROM_PAR to preserve each .par file EPHEM entry.")
    parser.add_argument("--seed", type=int, default=20260423)
    parser.add_argument("--output-suffix", default="", help="Optional suffix for output filenames.")
    parser.add_argument("--exclude-pulsars", nargs="*", default=[], help="Pulsar names to omit from this mechanical gate.")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    suffix = args.output_suffix.strip().replace("/", "_")
    output_tag = "baseline_powerlaw_pilot" if not suffix else f"baseline_powerlaw_pilot_{suffix}"
    result = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "task": "P2 timing-level baseline pilot",
        "status": "FAIL",
        "output_tag": output_tag,
        "arrays": args.arrays,
        "max_per_array": args.max_per_array,
        "components_gw": args.components_gw,
        "components_rn": args.components_rn,
        "include_red_noise": bool(args.include_red_noise),
        "ephem": args.ephem,
        "epta_variant": args.epta_variant,
        "seed": args.seed,
        "exclude_pulsars": sorted(set(args.exclude_pulsars)),
    }

    try:
        psrs, meta, skipped = load_pulsars(
            args.arrays,
            args.max_per_array,
            args.epta_variant,
            args.ephem,
            set(args.exclude_pulsars),
        )
        names = [m["name"] for m in meta]
        noise, provenance = load_noise_dicts(args.arrays, names, args.epta_variant)
        aliases = ppta_noise_aliases(psrs, meta, noise)
        pta, tspan = build_pta(psrs, meta, args.components_gw, args.components_rn, args.include_red_noise, noise)
        pta.set_default_params(noise)
        missing_white_noise = sorted(name for name in white_noise_params(list(pta.param_names)) if name not in noise)
        result.update(
            {
                "n_pulsars": len(psrs),
                "pulsars": meta,
                "skipped_pulsars": skipped,
                "n_noise_parameters": len(noise),
                "noise_provenance": provenance,
                "ppta_noise_aliases": aliases,
                "tspan_seconds": float(tspan),
                "tspan_years": float(tspan / 31557600.0),
                "ndim": len(pta.param_names),
                "param_names": list(pta.param_names),
                "missing_white_noise_parameters": missing_white_noise,
            }
        )
        if missing_white_noise:
            preview = ", ".join(missing_white_noise[:20])
            extra = "" if len(missing_white_noise) <= 20 else f", ... ({len(missing_white_noise)} total)"
            raise RuntimeError(f"Missing fixed white-noise parameters in supplied noise dictionaries: {preview}{extra}")
        x0 = make_initial_point(pta, noise, args.seed)
        lnprior = float(pta.get_lnprior(x0))
        lnlikelihood = float(pta.get_lnlikelihood(x0)) if math.isfinite(lnprior) else float("nan")
        result.update(
            {
                "status": "PASS" if math.isfinite(lnprior) and math.isfinite(lnlikelihood) else "FAIL",
                "lnprior": lnprior,
                "lnlikelihood": lnlikelihood,
            }
        )
    except Exception:
        result.update(
            {
                "error": traceback.format_exc(),
                "pulsars": result.get("pulsars", []),
                "skipped_pulsars": result.get("skipped_pulsars", []),
                "n_pulsars": result.get("n_pulsars", 0),
                "n_noise_parameters": result.get("n_noise_parameters", 0),
            }
        )

    json_path = OUT_DIR / f"{output_tag}.json"
    md_path = OUT_DIR / f"{output_tag}.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    write_markdown(md_path, result)
    print(json.dumps({"status": result["status"], "json": rel(json_path), "markdown": rel(md_path)}, indent=2))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
