#!/usr/bin/env python3
"""Read-only diagnostics for the remaining InPTA loader failures.

The script inspects staged public InPTA DR1 `.par` and `.tim` files and records
possible causes of the standard ENTERPRISE/libstempo loader failures.  It does
not edit timing models and does not construct replacement timing products.
"""

from __future__ import annotations

import csv
import json
import traceback
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPTA = ROOT / "data" / "InPTA_DR1" / "InPTA.DR1"
OUT = ROOT / "results" / "5pta_timing"


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def parse_fit_flag_anomalies(par: Path) -> list[dict]:
    rows = []
    for lineno, line in enumerate(par.read_text(errors="ignore").splitlines(), start=1):
        parts = line.split()
        if len(parts) < 3:
            continue
        key = parts[0]
        if not (key.startswith("DMXR") or key.startswith("DMXF")):
            continue
        if parts[-1] == "1":
            rows.append({"parfile": rel(par), "line": lineno, "parameter": key, "text": line.rstrip()})
    return rows


def tim_diagnostics(tim: Path) -> dict:
    lines = tim.read_text(errors="ignore").splitlines()
    first_nonempty = next((line.strip() for line in lines if line.strip()), "")
    data_like = []
    comment = 0
    malformed_c_prefix = []
    absolute_first_field = 0
    relative_parent_first_field = 0
    for lineno, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("C "):
            comment += 1
            continue
        if stripped.startswith("C") and not stripped.startswith("C "):
            malformed_c_prefix.append({"timfile": rel(tim), "line": lineno, "text": line.rstrip()})
        if stripped.startswith(("FORMAT", "MODE", "#")):
            continue
        data_like.append(stripped)
        first = stripped.split()[0]
        if first.startswith("/"):
            absolute_first_field += 1
        if first.startswith("../"):
            relative_parent_first_field += 1
    return {
        "timfile": rel(tim),
        "lines": len(lines),
        "first_nonempty": first_nonempty,
        "has_format_1": int(first_nonempty == "FORMAT 1" or any(line.strip() == "FORMAT 1" for line in lines[:3])),
        "has_mode_1": int(any(line.strip() == "MODE 1" for line in lines[:5])),
        "comment_rows": comment,
        "data_like_rows": len(data_like),
        "malformed_c_prefix_rows": len(malformed_c_prefix),
        "absolute_first_field_rows": absolute_first_field,
        "relative_parent_first_field_rows": relative_parent_first_field,
        "malformed_examples": malformed_c_prefix[:5],
    }


def try_enterprise_load(par: Path, tim: Path) -> dict:
    from enterprise.pulsar import Pulsar

    row = {"parfile": rel(par), "timfile": rel(tim)}
    for ephem in [None, "DE440", "DE436"]:
        key = "FROM_PAR" if ephem is None else ephem
        try:
            obj = Pulsar(str(par), str(tim), ephem=ephem)
            row[f"{key}_status"] = "PASS"
            row[f"{key}_n_toas"] = int(len(obj.toas))
            row[f"{key}_traceback_tail"] = ""
        except Exception:
            tb = traceback.format_exc()
            row[f"{key}_status"] = "FAIL"
            row[f"{key}_n_toas"] = ""
            row[f"{key}_traceback_tail"] = "\n".join(tb.splitlines()[-8:])
    return row


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    targets = []
    for psr_dir in sorted(p for p in INPTA.glob("*") if p.is_dir()):
        psr = psr_dir.name
        for band in ["NB", "WB"]:
            par = psr_dir / f"{psr}.{band}.par"
            tim = psr_dir / f"{psr}.{band}.tim"
            if par.exists() and tim.exists():
                targets.append((psr, band, par, tim))

    par_rows = []
    tim_rows = []
    load_rows = []
    malformed_rows = []
    for psr, band, par, tim in targets:
        for row in parse_fit_flag_anomalies(par):
            row.update({"pulsar": psr, "band": band})
            par_rows.append(row)
        trow = tim_diagnostics(tim)
        trow.update({"pulsar": psr, "band": band})
        tim_rows.append(trow)
        malformed_rows.extend(dict(ex, pulsar=psr, band=band) for ex in trow.pop("malformed_examples", []))

    focused = [
        ("J0751+1807", "NB", INPTA / "J0751+1807" / "J0751+1807.NB.par", INPTA / "J0751+1807" / "J0751+1807.NB.tim"),
        ("J0751+1807", "WB", INPTA / "J0751+1807" / "J0751+1807.WB.par", INPTA / "J0751+1807" / "J0751+1807.WB.tim"),
        ("J0437-4715", "NB", INPTA / "J0437-4715" / "J0437-4715.NB.par", INPTA / "J0437-4715" / "J0437-4715.NB.tim"),
        ("J0437-4715", "WB", INPTA / "J0437-4715" / "J0437-4715.WB.par", INPTA / "J0437-4715" / "J0437-4715.WB.tim"),
    ]
    for psr, band, par, tim in focused:
        row = try_enterprise_load(par, tim)
        row.update({"pulsar": psr, "band": band})
        load_rows.append(row)

    write_csv(OUT / "inpta_par_fitflag_diagnostics.csv", par_rows, ["pulsar", "band", "parfile", "line", "parameter", "text"])
    write_csv(
        OUT / "inpta_tim_format_diagnostics.csv",
        tim_rows,
        [
            "pulsar",
            "band",
            "timfile",
            "lines",
            "first_nonempty",
            "has_format_1",
            "has_mode_1",
            "comment_rows",
            "data_like_rows",
            "malformed_c_prefix_rows",
            "absolute_first_field_rows",
            "relative_parent_first_field_rows",
        ],
    )
    write_csv(OUT / "inpta_tim_malformed_c_prefix.csv", malformed_rows, ["pulsar", "band", "timfile", "line", "text"])
    write_csv(
        OUT / "inpta_focused_enterprise_loads.csv",
        load_rows,
        [
            "pulsar",
            "band",
            "parfile",
            "timfile",
            "FROM_PAR_status",
            "FROM_PAR_n_toas",
            "FROM_PAR_traceback_tail",
            "DE440_status",
            "DE440_n_toas",
            "DE440_traceback_tail",
            "DE436_status",
            "DE436_n_toas",
            "DE436_traceback_tail",
        ],
    )

    nb_fitflags = [row for row in par_rows if row["band"] == "NB"]
    wb_no_format = [row for row in tim_rows if row["band"] == "WB" and not row["has_format_1"]]
    lines = [
        "# InPTA Loader Failure Diagnostics",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Scope",
        "",
        "- Read-only inspection of staged public InPTA DR1 par/tim files.",
        "- No `.par` or `.tim` file is edited.",
        "- Diagnostics identify likely loader issues; they do not define a production timing policy.",
        "",
        "## Summary",
        "",
        f"- InPTA par/tim targets inspected: `{len(targets)}`.",
        f"- DMXR/DMXF fit-flag anomalies found: `{len(par_rows)}`.",
        f"- NB fit-flag anomalies found: `{len(nb_fitflags)}`.",
        f"- WB tim files without `FORMAT 1`: `{len(wb_no_format)}`.",
        f"- Malformed `C`-prefix tim rows: `{len(malformed_rows)}`.",
        "",
        "## Focused Load Attempts",
        "",
        "| pulsar | band | FROM_PAR | DE440 | DE436 |",
        "|---|---|---|---|---|",
    ]
    for row in load_rows:
        lines.append(
            f"| `{row['pulsar']}` | `{row['band']}` | `{row['FROM_PAR_status']}` | `{row['DE440_status']}` | `{row['DE436_status']}` |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- `J0751+1807.NB.par` contains a fitted `DMXR1_0007` line; libstempo reports `param_dmxr1` during the failing full loader run, so this is the leading read-only explanation for the NB design-matrix failure.",
            "- The WB `.tim` files lack the `FORMAT 1` header used by the NB files and all focused WB ENTERPRISE load attempts fail before usable TOAs are exposed.",
            "- A production fix should either use an official InPTA loading/preprocessing recipe or define a documented exclusion/preprocessing policy; this diagnostic intentionally does not patch the timing files.",
            "",
            "## Output Files",
            "",
            "- `results/5pta_timing/inpta_par_fitflag_diagnostics.csv`",
            "- `results/5pta_timing/inpta_tim_format_diagnostics.csv`",
            "- `results/5pta_timing/inpta_tim_malformed_c_prefix.csv`",
            "- `results/5pta_timing/inpta_focused_enterprise_loads.csv`",
            "",
        ]
    )
    (OUT / "inpta_loader_failure_diagnostic.md").write_text("\n".join(lines))
    (OUT / "inpta_loader_failure_diagnostic.json").write_text(
        json.dumps(
            {
                "generated": datetime.now().isoformat(timespec="seconds"),
                "targets": len(targets),
                "fitflag_anomalies": len(par_rows),
                "nb_fitflag_anomalies": len(nb_fitflags),
                "wb_without_format_1": len(wb_no_format),
                "malformed_c_prefix_rows": len(malformed_rows),
                "report": rel(OUT / "inpta_loader_failure_diagnostic.md"),
            },
            indent=2,
        )
        + "\n"
    )
    print({"targets": len(targets), "fitflag_anomalies": len(par_rows), "wb_without_format_1": len(wb_no_format)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
