#!/usr/bin/env python3
"""Static gate for the local PRL submission package.

This script does not run TeX.  It records whether the package is structurally
ready for the final compile/arbiter pass and makes unresolved human blockers
explicit.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "results/T2_NG15yr/prl_package_static_gate.json"
OUT_MD = ROOT / "results/T2_NG15yr/prl_package_static_gate.md"

TITLE = "Calibrated PTA Evidence Favors Curved Spectra but Not a Unique Nanohertz Source"

REQUIRED_FILES = {
    "revtex": "theory/paper_prl_submission.tex",
    "compressed_md": "theory/paper_prl_compressed_draft.md",
    "supplement_md": "theory/prl_supplement_draft.md",
    "supplement_tex": "theory/prl_supplement.tex",
    "package_manifest": "theory/PRL_submission_package_manifest.md",
    "tex_static_check": "theory/T3.11_PRL_tex_static_check.md",
    "formal_compile_report": "theory/T3.12_PRL_formal_compile_report.md",
    "human_gate": "theory/PRL_submission_human_gate.md",
    "cover_letter": "theory/PRL_cover_letter_draft.md",
    "figure": "results/T2_NG15yr/figures/prl_decisive_evidence_figure.pdf",
    "bridge_figure": "results/T2_NG15yr/figures/prl_bridge_evidence_figure.pdf",
    "bridge_summary": "results/prl_reference_bridge/P0_P6_run_summary.md",
    "revtex_pdf": "theory/pdf/revtex/paper_prl_submission.pdf",
    "supplement_pdf": "theory/pdf/revtex/prl_supplement.pdf",
    "reproducibility": "REPRODUCIBILITY.md",
    "environment": "environment.yml",
}

JSON_RESULTS = {
    "official_density_robustness": "results/T2_NG15yr/bayes_factors/prl_official_density_robustness.json",
    "H4_env_robustness": "results/T2_NG15yr/bayes_factors/prl_H4_env_robustness.json",
    "H4_H5_H6": "results/T2_NG15yr/bayes_factors/prl_H4_H5_H6_experiments.json",
    "H7_astro_family": "results/T2_NG15yr/bayes_factors/prl_H7_astrophysical_family.json",
    "H9_bin_driver": "results/T2_NG15yr/bayes_factors/prl_H9_bin_driver_analysis.json",
    "H10_systematic_envelope": "results/T2_NG15yr/bayes_factors/prl_H10_systematic_envelope.json",
    "ti_qmc_crosscheck": "results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json",
    "car_null_calibration": "results/T2_NG15yr/covariance/car_null_calibration.json",
    "rhat": "results/T2_NG15yr/T2_multichain_Rhat.json",
}

TEX_ENGINES = ["pdflatex", "latexmk", "tectonic", "xelatex", "lualatex"]
HARD_PLACEHOLDERS = ["[CODE-RELEASE-URL-TBD]"]
SOFT_PLACEHOLDERS = ["[Acknowledgments, grants, and computing-resource text TBD.]"]


def word_count(text: str) -> int:
    return len(text.split())


def cite_report(tex: str) -> dict:
    cites = set()
    for body in re.findall(r"\\cite\{([^}]*)\}", tex):
        cites.update(key.strip() for key in body.split(",") if key.strip())
    bibs = set(re.findall(r"\\bibitem\{([^}]*)\}", tex))
    return {
        "cite_count": len(cites),
        "bibitem_count": len(bibs),
        "missing_bibitems": sorted(cites - bibs),
        "unused_bibitems": sorted(bibs - cites),
    }


def env_report(tex: str) -> dict:
    stack: list[str] = []
    for match in re.finditer(r"\\(begin|end)\{([^}]*)\}", tex):
        typ, name = match.groups()
        if typ == "begin":
            stack.append(name)
        elif not stack or stack[-1] != name:
            return {
                "balanced": False,
                "mismatch": {"type": typ, "name": name, "stack_tail": stack[-5:]},
            }
        else:
            stack.pop()
    return {"balanced": not stack, "remaining_stack": stack}


def process_scan() -> dict:
    pattern = re.compile(r"prl_H4|prl_H5|prl_H6|NG15yr|bayes_factors|PTArcade|PTMCMC|dynesty")
    python_process = re.compile(r"python|Python", re.IGNORECASE)
    proc = subprocess.run(
        ["ps", "-axo", "pid,stat,etime,command"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    lines = []
    for line in proc.stdout.splitlines():
        if not pattern.search(line):
            continue
        if not python_process.search(line):
            continue
        if "prl_package_static_gate.py" in line:
            continue
        if "-m py_compile" in line or "python3 - <<" in line:
            continue
        lines.append(line.strip())
    return {"matching_processes": lines, "count": len(lines)}


def main() -> None:
    generated = datetime.now().isoformat(timespec="seconds")
    file_status = {
        label: {"path": path, "exists": (ROOT / path).exists()}
        for label, path in REQUIRED_FILES.items()
    }

    tex_path = ROOT / REQUIRED_FILES["revtex"]
    tex = tex_path.read_text() if tex_path.exists() else ""
    compressed_path = ROOT / REQUIRED_FILES["compressed_md"]
    compressed = compressed_path.read_text() if compressed_path.exists() else ""
    supplement_path = ROOT / REQUIRED_FILES["supplement_md"]
    supplement = supplement_path.read_text() if supplement_path.exists() else ""
    supplement_tex_path = ROOT / REQUIRED_FILES["supplement_tex"]
    supplement_tex = supplement_tex_path.read_text() if supplement_tex_path.exists() else ""
    cover_path = ROOT / REQUIRED_FILES["cover_letter"]
    cover = cover_path.read_text() if cover_path.exists() else ""

    json_status = {}
    for label, relpath in JSON_RESULTS.items():
        path = ROOT / relpath
        status = {"path": relpath, "exists": path.exists(), "complete": None, "readable": False}
        if path.exists():
            try:
                data = json.loads(path.read_text())
                status["readable"] = True
                status["complete"] = data.get("complete", True)
            except Exception as exc:  # noqa: BLE001 - record failure in gate output
                status["error"] = repr(exc)
        json_status[label] = status

    placeholders = {
        "hard": [p for p in HARD_PLACEHOLDERS if p in tex or p in compressed or p in supplement_tex or p in cover],
        "soft": [p for p in SOFT_PLACEHOLDERS if p in tex or p in compressed or p in supplement_tex],
    }

    main_citations = cite_report(tex)
    supplement_citations = cite_report(supplement_tex)
    main_env = env_report(tex)
    supplement_env = env_report(supplement_tex)
    tex_engines = {engine: shutil.which(engine) for engine in TEX_ENGINES}
    checks = {
        "title_present_tex": TITLE in tex,
        "title_present_compressed": TITLE in compressed,
        "title_present_supplement_tex": TITLE in supplement_tex,
        "all_required_files_exist": all(item["exists"] for item in file_status.values()),
        "all_json_readable": all(item["readable"] for item in json_status.values()),
        "all_json_complete": all(item["complete"] is True for item in json_status.values()),
        "no_missing_bibitems": not main_citations["missing_bibitems"],
        "no_unused_bibitems": not main_citations["unused_bibitems"],
        "supplement_no_missing_bibitems": not supplement_citations["missing_bibitems"],
        "supplement_no_unused_bibitems": not supplement_citations["unused_bibitems"],
        "latex_env_balanced": main_env["balanced"],
        "supplement_latex_env_balanced": supplement_env["balanced"],
        "no_hard_placeholders": not placeholders["hard"],
        "no_soft_placeholders": not placeholders["soft"],
        "tex_engine_available": any(tex_engines.values()),
        "no_matching_background_processes": process_scan()["count"] == 0,
    }
    structural_pass = all(
        checks[key]
        for key in [
            "title_present_tex",
            "title_present_compressed",
            "title_present_supplement_tex",
            "all_required_files_exist",
            "all_json_readable",
            "all_json_complete",
            "no_missing_bibitems",
            "no_unused_bibitems",
            "supplement_no_missing_bibitems",
            "supplement_no_unused_bibitems",
            "latex_env_balanced",
            "supplement_latex_env_balanced",
            "no_matching_background_processes",
        ]
    )
    submission_ready = (
        structural_pass
        and checks["no_hard_placeholders"]
        and checks["no_soft_placeholders"]
        and checks["tex_engine_available"]
    )

    report = {
        "generated": generated,
        "title": TITLE,
        "checks": checks,
        "structural_pass": structural_pass,
        "submission_ready": submission_ready,
        "file_status": file_status,
        "json_status": json_status,
        "word_counts": {
            "revtex_source": word_count(tex),
            "compressed_markdown": word_count(compressed),
            "supplement_markdown": word_count(supplement),
            "supplement_tex_source": word_count(supplement_tex),
        },
        "citations": main_citations,
        "supplement_citations": supplement_citations,
        "latex_environments": main_env,
        "supplement_latex_environments": supplement_env,
        "placeholders": placeholders,
        "tex_engines": tex_engines,
        "process_scan": process_scan(),
        "blockers": [
            "public code-release URL is unresolved" if placeholders["hard"] else None,
            "soft placeholders remain" if placeholders["soft"] else None,
            "no TeX engine available" if not checks["tex_engine_available"] else None,
        ],
    }
    report["blockers"] = [item for item in report["blockers"] if item]

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    lines = [
        "# PRL Package Static Gate",
        "",
        f"**Generated**: {generated}",
        f"**Structural pass**: `{structural_pass}`",
        f"**Submission ready**: `{submission_ready}`",
        "",
        "## Checks",
        "",
        "| Check | Pass |",
        "|---|---:|",
    ]
    for key, value in checks.items():
        lines.append(f"| `{key}` | `{value}` |")
    lines.extend(
        [
            "",
            "## Word Counts",
            "",
            "| File | words |",
            "|---|---:|",
            f"| REVTeX source | `{report['word_counts']['revtex_source']}` |",
            f"| compressed Markdown | `{report['word_counts']['compressed_markdown']}` |",
            f"| supplement Markdown | `{report['word_counts']['supplement_markdown']}` |",
            f"| supplement TeX source | `{report['word_counts']['supplement_tex_source']}` |",
            "",
            "## Blockers",
            "",
        ]
    )
    if report["blockers"]:
        lines.extend(f"- {item}" for item in report["blockers"])
    else:
        lines.append("- none")
    lines.extend(
        [
            "",
            "## TeX Engines",
            "",
            "| Engine | path |",
            "|---|---|",
        ]
    )
    for engine, path in tex_engines.items():
        lines.append(f"| `{engine}` | `{path or ''}` |")
    lines.extend(
        [
            "",
            "## Background Process Scan",
            "",
            f"Matching process count: `{report['process_scan']['count']}`",
        ]
    )
    for line in report["process_scan"]["matching_processes"]:
        lines.append(f"- `{line}`")
    OUT_MD.write_text("\n".join(lines) + "\n")
    print(f"saved: {OUT_JSON}")
    print(f"saved: {OUT_MD}")
    print(f"structural_pass={structural_pass} submission_ready={submission_ready}")


if __name__ == "__main__":
    main()
