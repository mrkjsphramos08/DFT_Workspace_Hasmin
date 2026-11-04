#!/usr/bin/env python3
"""Summarize VASP DFT relaxation logs with Python 3.8+; no extra packages.

Examples:
  python3 vasp_convergence.py /path/to/run
  python3 vasp_convergence.py . --save convergence
  python3 vasp_convergence.py OUTCAR --poscar POSCAR --save convergence
  python3 vasp_convergence.py vasp.out --incar INCAR

Reads one OUTCAR, OSZICAR, or stdout snapshot. Never edits VASP inputs and
never reads POTCAR. Full free/fixed atom masks are supported. Partial
directional constraints are deliberately reported as unsupported.

Outputs include plain text, JSON, CSV (ionic and SCF), and optionally
Markdown (.md) for sharing results.
"""
import argparse
import csv
import gzip
import io
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

VERSION = "1.0.0"
NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
KEYS = {"SYSTEM", "EDIFF", "EDIFFG", "NELM", "NSW", "IBRION", "ISIF",
        "ISPIN", "ENCUT", "SIGMA", "ISMEAR", "POTIM", "LDIPOL", "IDIPOL",
        "DIPOL", "LORBIT", "NUPDOWN", "IVDW", "POMASS", "ICHARG", "ISTART"}
FLOAT_KEYS = {"EDIFF", "EDIFFG", "ENCUT", "SIGMA", "POTIM", "NUPDOWN"}
INT_KEYS = {"NELM", "NSW", "IBRION", "ISIF", "ISPIN", "ISMEAR", "IDIPOL",
            "LORBIT", "IVDW", "ICHARG", "ISTART"}


def number(value):
    result = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError("Non-finite number")
    return result


def numbers(line):
    return [number(x) for x in re.findall(NUM, line)]


def read_text(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def settings_from_line(line):
    result = {}
    line = re.split(r"[!#]", line, maxsplit=1)[0]
    for match in re.finditer(r"\b([A-Z][A-Z0-9_]*)\s*=\s*([^;]*?)(?=\s+[A-Z][A-Z0-9_]*\s*=|;|$)", line):
        key, raw = match.group(1), match.group(2).strip()
        if key not in KEYS:
            continue
        try:
            if key in FLOAT_KEYS:
                result[key] = number(re.match(NUM, raw).group())
            elif key in INT_KEYS:
                result[key] = int(number(re.match(NUM, raw).group()))
            elif key == "LDIPOL":
                result[key] = raw.upper().strip(".").startswith("T")
            elif key in {"DIPOL", "POMASS"}:
                # Echoed OUTCAR POMASS entries are separate from POTCAR metadata.
                result[key] = numbers(raw)
            else:
                result[key] = raw
        except (ValueError, AttributeError):
            pass
    return result


def read_settings(path):
    result = {}
    if path:
        for line in read_text(path).splitlines():
            result.update(settings_from_line(line))
    return result


def read_poscar(path):
    """Read counts/masks only; preserve arbitrary species labels and atom order."""
    lines = read_text(path).splitlines()
    if len(lines) < 8:
        raise ValueError("POSCAR is too short")
    cursor = 5
    if not all(re.fullmatch(r"\d+", x) for x in lines[cursor].split()):
        cursor += 1  # VASP 5 species line
    counts = [int(x) for x in lines[cursor].split()]
    if not counts or any(x < 0 for x in counts) or sum(counts) == 0:
        raise ValueError("Invalid POSCAR atom counts")
    cursor += 1
    selective = lines[cursor].strip().lower().startswith("s")
    if selective:
        cursor += 1
    mode = lines[cursor].strip().lower()
    if not mode or mode[0] not in "dck":
        raise ValueError("POSCAR coordinate mode is not recognized")
    cursor += 1
    flags = []
    for row in lines[cursor:cursor + sum(counts)]:
        fields = row.split()
        if len(fields) < (6 if selective else 3):
            raise ValueError("Incomplete POSCAR coordinate/flag row")
        [number(x) for x in fields[:3]]
        atom_flags = [True, True, True]
        if selective:
            if any(x.upper() not in {"T", "F"} for x in fields[3:6]):
                raise ValueError("Invalid POSCAR selective-dynamics flags")
            atom_flags = [x.upper() == "T" for x in fields[3:6]]
        flags.append(atom_flags)
    if len(flags) != sum(counts):
        raise ValueError("POSCAR atom rows do not match counts")
    return {"counts": counts, "atoms": sum(counts), "selective": selective,
            "flags": flags, "partial": any(any(f) and not all(f) for f in flags),
            "free_atoms": sum(all(f) for f in flags),
            "fixed_atoms": sum(not any(f) for f in flags)}


def new_step(index):
    return {"evaluation": index, "ionic_complete": False, "scf": [],
            "scf_convergence": "unconfirmed", "final_F_eV": None,
            "E0_eV": None, "energy_without_entropy_eV": None,
            "magnetization_muB": None, "raw_max_force_eV_A": None,
            "free_max_force_eV_A": None, "raw_max_force_atom": None,
            "free_max_force_atom": None, "force_drift_eV_A": None,
            "force_drift_norm_eV_A": None, "pressure_kB": None,
            "dispersion_eV": None, "dipole_e_A": None,
            "force_table_complete": False, "notes": []}


def add_warning(warnings, line_number, message):
    warnings.append({"line": line_number, "message": message.strip()})


def observe_messages(line, lineno, meta):
    lower = line.lower()
    if "warning" in lower or "zbrent:" in lower or "brmix:" in lower:
        add_warning(meta["warnings"], lineno, line)
    if re.search(r"VERY BAD NEWS|internal error|MPI_ABORT|segmentation fault|SIGSEGV|"
                 r"ZBRENT: fatal|Call to ZHEGV failed|BRMIX: very serious", line, re.I):
        meta["errors"].append({"line": lineno, "message": line.strip()})
    if "reached required accuracy - stopping structural energy minimisation" in lower:
        meta["structural_accuracy_message"] = True
    if "general timing and accounting" in lower:
        meta["normal_footer_seen"] = True
    if "elapsed time (sec):" in lower:
        vals = numbers(line)
        if vals:
            meta["elapsed_seconds"] = vals[-1]
    if re.search(r"LSTOP\s*=\s*\.?T|LABORT\s*=\s*\.?T|soft stop encountered|STOPCAR found", line, re.I):
        meta["stop_request_seen"] = True


def force_stats(step, rows, mask):
    if not rows:
        return
    norms = [math.sqrt(sum(v * v for v in row[3:6])) for row in rows]
    step["raw_max_force_eV_A"] = max(norms)
    step["raw_max_force_atom"] = norms.index(max(norms)) + 1
    step["force_table_complete"] = True
    if mask and len(mask["flags"]) == len(rows) and not mask["partial"]:
        free = [(i + 1, norm) for i, (norm, f) in enumerate(zip(norms, mask["flags"])) if all(f)]
        if free:
            atom, value = max(free, key=lambda item: item[1])
            step["free_max_force_eV_A"] = value
            step["free_max_force_atom"] = atom
        else:
            step["notes"].append("No free atoms; no meaningful free-atom relaxation test")


def parse_outcar(lines, meta, mask):
    steps = []
    current = None
    scf = None
    seen = set()
    header_settings = {}
    fallback_settings = {}
    in_incar = False
    final_block = False
    nions = None
    force_rows = None
    force_bad = False
    dataset_masses = []
    for lineno, line in enumerate(lines, 1):
        observe_messages(line, lineno, meta)
        if line.strip() == "INCAR:":
            in_incar = True
            continue
        if in_incar and line.strip().startswith("POTCAR:"):
            in_incar = False
        if in_incar:
            header_settings.update(settings_from_line(line))
        if not steps:
            if "ZVAL" in line and "POMASS" in line:
                m = re.search(r"POMASS\s*=\s*(" + NUM + ")", line)
                if m:
                    dataset_masses.append(number(m.group(1)))
            elif not in_incar:
                for key, value in settings_from_line(line).items():
                    fallback_settings.setdefault(key, value)
            m = re.search(r"NIONS\s*=\s*(\d+)", line)
            if m:
                nions = int(m.group(1))
            if "ions per type" in line:
                meta["ions_per_type"] = [int(x) for x in line.split("=", 1)[1].split()]
                nions = sum(meta["ions_per_type"])
            if "using selective dynamics" in line.lower():
                meta["selective_dynamics_in_log"] = True
        m = re.search(r"\bIteration\s+(\d+)\(\s*(\d+)\s*\)", line)
        if m:
            index, iteration = map(int, m.groups())
            if (index, iteration) in seen:
                raise ValueError("Repeated OUTCAR iteration labels: split concatenated runs before analysis")
            seen.add((index, iteration))
            if current is None or current["evaluation"] != index:
                if current is not None and index <= current["evaluation"]:
                    raise ValueError("OUTCAR ionic numbering resets: analyze separate runs separately")
                current = new_step(index)
                steps.append(current)
            scf = {"evaluation": index, "iteration": iteration,
                   "energy_eV": None, "dE_eV": None, "d_eps_eV": None,
                   "rms_charge": None, "magnetization_muB": None,
                   "loop_seconds": None, "line": lineno}
            current["scf"].append(scf)
            final_block = False
            force_rows = None
            continue
        if current is None:
            continue
        if "aborting loop because EDIFF is reached" in line:
            current["scf_convergence"] = "reached (OUTCAR message)"
        if scf is not None:
            m = re.search(r"total energy-change \(2\. order\)\s*:\s*(" + NUM + r")\s*\(\s*(" + NUM + r")", line)
            if m:
                scf["dE_eV"], scf["d_eps_eV"] = map(number, m.groups())
            m = re.search(r"number of electron.*?magnetization\s+(" + NUM + ")", line)
            if m:
                scf["magnetization_muB"] = number(m.group(1))
                current["magnetization_muB"] = scf["magnetization_muB"]
            m = re.search(r"rms\(total\)\s*=\s*(" + NUM + ")", line)
            if m:
                scf["rms_charge"] = number(m.group(1))
            m = re.search(r"\bLOOP:\s+cpu time.*?real time\s+(" + NUM + ")", line)
            if m:
                scf["loop_seconds"] = number(m.group(1))
        if "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in line:
            final_block = True
            continue
        m = re.search(r"free\s+energy\s+TOTEN\s*=\s*(" + NUM + ")", line)
        if m:
            if final_block:
                current["final_F_eV"] = number(m.group(1))
                current["ionic_complete"] = True
            elif scf is not None:
                scf["energy_eV"] = number(m.group(1))
        if final_block and "energy(sigma->0)" in line:
            m = re.search(r"energy\s+without entropy\s*=\s*(" + NUM + r").*?energy\(sigma->0\)\s*=\s*(" + NUM + ")", line)
            if m:
                current["energy_without_entropy_eV"], current["E0_eV"] = map(number, m.groups())
            final_block = False
        m = re.search(r"Edisp \(eV\)\s+(" + NUM + ")", line)
        if m:
            current["dispersion_eV"] = number(m.group(1))
        m = re.search(r"external pressure\s*=\s*(" + NUM + ")", line)
        if m:
            current["pressure_kB"] = number(m.group(1))
        if line.strip().startswith("dipolmoment"):
            vals = numbers(line)
            if len(vals) == 3:
                current["dipole_e_A"] = vals
        if "POSITION" in line and "TOTAL-FORCE" in line:
            force_rows, force_bad = [], False
            continue
        if force_rows is not None:
            stripped = line.strip()
            if not stripped or set(stripped) == {"-"}:
                if force_rows:
                    if (nions is None or len(force_rows) != nions) or force_bad:
                        current["notes"].append("Incomplete or malformed force table; no force result used")
                    else:
                        force_stats(current, force_rows, mask)
                    force_rows = None
                continue
            fields = line.split()
            try:
                row = [number(x) for x in fields[:6]]
                if len(row) != 6:
                    raise ValueError("short row")
                force_rows.append(row)
            except ValueError:
                force_bad = True
                force_rows = None
                current["notes"].append("Malformed force table; no force result used")
        if "total drift:" in line:
            vals = numbers(line.split("total drift:", 1)[1])
            if len(vals) == 3:
                current["force_drift_eV_A"] = vals
                current["force_drift_norm_eV_A"] = math.sqrt(sum(x*x for x in vals))
    if force_rows is not None:
        current["notes"].append("File ends inside force table; no force result used")
    meta["atom_count"] = nions
    settings = dict(fallback_settings)
    settings.update(header_settings)  # preserve exact INCAR values rather than rounded echo
    meta["settings"] = settings
    if dataset_masses and "POMASS" in header_settings:
        requested = header_settings["POMASS"]
        if len(requested) == len(dataset_masses):
            mismatches = [{"species_index": i + 1, "INCAR_mass": a, "potential_mass": b}
                          for i, (a, b) in enumerate(zip(requested, dataset_masses)) if abs(a-b) > 1e-5]
            meta["mass_overrides"] = mismatches
    return steps


def parse_compact(lines, meta, settings):
    """OSZICAR and stdout share DAV/RMM/CG electronic and F= ionic rows."""
    steps, current = [], new_step(1)
    last_index = 0
    for lineno, line in enumerate(lines, 1):
        observe_messages(line, lineno, meta)
        m = re.match(r"\s*([A-Z][A-Z0-9_-]*)\s*:\s*(\d+)\s+(.*)", line)
        if m and m.group(1) in {"DAV", "RMM", "CG", "DMP", "RMM-DIIS"}:
            fields = m.group(3).split()
            if len(fields) < 3:
                current["notes"].append("Incomplete electronic row at line %d" % lineno)
                continue
            try:
                row = {"evaluation": current["evaluation"], "iteration": int(m.group(2)),
                       "energy_eV": number(fields[0]), "dE_eV": number(fields[1]),
                       "d_eps_eV": number(fields[2]), "rms_charge": None,
                       "magnetization_muB": None, "loop_seconds": None, "line": lineno}
                if len(fields) >= 6:
                    row["rms_charge"] = number(fields[5])
                if current["scf"] and row["iteration"] <= current["scf"][-1]["iteration"]:
                    raise ValueError("Electronic numbering reset before an ionic summary")
                current["scf"].append(row)
            except ValueError as exc:
                if "numbering reset" in str(exc):
                    raise
                current["notes"].append("Malformed electronic row at line %d" % lineno)
        m = re.match(r"\s*(\d+)\s+F=\s*(" + NUM + r")\s+E0=\s*(" + NUM + ")", line)
        if m:
            index = int(m.group(1))
            if index <= last_index:
                raise ValueError("Ionic numbering resets: split concatenated logs before analysis")
            last_index = index
            current["evaluation"] = index
            for row in current["scf"]:
                row["evaluation"] = index
            current["ionic_complete"] = True
            current["final_F_eV"], current["E0_eV"] = number(m.group(2)), number(m.group(3))
            mag = re.search(r"mag=\s*(" + NUM + ")", line)
            if mag:
                current["magnetization_muB"] = number(mag.group(1))
            steps.append(current)
            current = new_step(index + 1)
    if current["scf"] or current["notes"]:
        steps.append(current)
    meta["settings"] = settings
    meta["limitations"].append("Compact log: forces, pressure, and OUTCAR completion footer are unavailable")
    return steps


def finish_analysis(steps, meta, args, mask):
    settings = meta["settings"]
    ediff = settings.get("EDIFF")
    ediffg = settings.get("EDIFFG")
    nelm = settings.get("NELM")
    previous = None
    flags = []
    for step in steps:
        recorded = [row for row in step["scf"] if row["energy_eV"] is not None]
        step["scf_iterations_started"] = max((r["iteration"] for r in step["scf"]), default=0)
        step["scf_iterations_with_energy"] = len(recorded)
        last = recorded[-1] if recorded else None
        step["last_scf_energy_eV"] = last["energy_eV"] if last else None
        step["last_scf_dE_eV"] = last["dE_eV"] if last else None
        step["last_scf_d_eps_eV"] = last["d_eps_eV"] if last else None
        step["last_scf_rms_charge"] = last["rms_charge"] if last else None
        step["energy_criterion_satisfied"] = None
        if ediff is not None and ediff > 0 and last and last["dE_eV"] is not None and last["d_eps_eV"] is not None:
            passed = abs(last["dE_eV"]) < ediff and abs(last["d_eps_eV"]) < ediff
            step["energy_criterion_satisfied"] = passed
            if step["scf_convergence"] == "unconfirmed" and step["ionic_complete"]:
                step["scf_convergence"] = "inferred from energy criterion" if passed else "not demonstrated by energy criterion"
        if nelm and step["scf_iterations_started"] >= nelm and not step["scf_convergence"].startswith("reached"):
            step["notes"].append("NELM reached without explicit convergence confirmation")
            flags.append("Evaluation %d reached NELM without explicit convergence confirmation" % step["evaluation"])
        step["delta_F_eV"] = None
        step["free_force_below_target"] = None
        if ediffg is not None and ediffg < 0 and step["free_max_force_eV_A"] is not None:
            step["free_force_below_target"] = step["free_max_force_eV_A"] < abs(ediffg)
        if step["ionic_complete"] and previous is not None:
            delta = step["final_F_eV"] - previous["final_F_eV"]
            step["delta_F_eV"] = delta
            if delta > args.energy_jump:
                flags.append("Evaluation %d: F increased by %.6f eV (review threshold %.3f eV/cell)" % (step["evaluation"], delta, args.energy_jump))
            now, before = step["free_max_force_eV_A"], previous["free_max_force_eV_A"]
            if now is not None and before is not None and now > args.force_warning and now > 3 * before:
                flags.append("Evaluation %d: free-atom force jumped from %.6f to %.6f eV/A" % (step["evaluation"], before, now))
        if step["force_drift_norm_eV_A"] is not None and step["force_drift_norm_eV_A"] > args.drift_warning:
            flags.append("Evaluation %d: reported force-drift norm %.6f eV/A exceeds review threshold %.3f" % (step["evaluation"], step["force_drift_norm_eV_A"], args.drift_warning))
        if step["ionic_complete"]:
            previous = step
    meta["review_flags"] = flags
    meta["heuristic_thresholds"] = {"positive_delta_F_eV_per_cell": args.energy_jump,
                                    "force_spike_minimum_eV_A": args.force_warning,
                                    "force_spike_factor": 3,
                                    "force_drift_eV_A": args.drift_warning}
    complete = [s for s in steps if s["ionic_complete"]]
    latest = complete[-1] if complete else None
    meta["completed_ionic_evaluations"] = len(complete)
    meta["electronic_cycles_confirmed_by_OUTCAR"] = sum(s["scf_convergence"].startswith("reached") for s in steps)
    meta["last_completed_evaluation"] = latest["evaluation"] if latest else None
    meta["latest_observed_evaluation"] = steps[-1]["evaluation"] if steps else None
    meta["unfinished_evaluation_present"] = bool(steps and not steps[-1]["ionic_complete"])
    static = settings.get("NSW") == 0 or settings.get("IBRION") == -1
    meta["calculation_mode"] = "static" if static else "molecular dynamics" if settings.get("IBRION") == 0 else "relaxation or unknown"
    if static:
        meta["ionic_status"] = "Not applicable: static calculation; forces do not imply relaxation"
    elif meta["structural_accuracy_message"]:
        meta["ionic_status"] = "VASP reports structural convergence"
        if latest and latest["free_force_below_target"] is False:
            meta["ionic_status"] += "; supplied-mask force check conflicts: investigate"
        if settings.get("ISIF", 2) != 2:
            meta["limitations"].append("Cell degrees of freedom are active or different: atomic forces alone do not independently verify cell convergence")
    elif meta["unfinished_evaluation_present"]:
        meta["ionic_status"] = "Not confirmed: log ends during an unfinished evaluation"
    elif settings.get("NSW", 0) > 0 and len(complete) >= settings["NSW"]:
        meta["ionic_status"] = "Not confirmed: ionic step limit reached without structural-accuracy message"
    else:
        meta["ionic_status"] = "Not confirmed: no structural-accuracy message found"
    if meta["errors"]:
        meta["ending_status"] = "Error text detected; inspect the reported lines"
    elif meta["normal_footer_seen"]:
        meta["ending_status"] = "Normal accounting footer observed; this does not establish ionic convergence"
    else:
        meta["ending_status"] = "No normal accounting footer observed: running, stopped, or truncated cannot be distinguished from this snapshot"
    if not mask:
        meta["limitations"].append("No usable POSCAR mask: raw force maxima are NOT automatically called free-atom force maxima")
    elif mask["partial"]:
        meta["limitations"].append("Partial directional T/F constraints: free-force projection is not implemented; raw forces only")
    if meta["source_format"] != "outcar":
        meta["limitations"].append("Electronic convergence from energy columns is inferred, not an explicit OUTCAR confirmation")
    if not settings.get("EDIFF"):
        meta["limitations"].append("Positive EDIFF unavailable: energy-column convergence cannot be evaluated")
    meta["limitations"].append("A compatible POSCAR atom order is required; matching counts alone cannot prove order provenance")
    meta["limitations"].append("This is log analysis, not validation of the surface model, magnetic ground state, or convergence of cutoff/k-points/vacuum")


def fmt(value, precision=6):
    return "--" if value is None else ("%.*f" % (precision, value))


def report_markdown(steps, meta, last_count):
    """Generate a Markdown formatted report."""
    settings = meta["settings"]
    lines = [
        f"# VASP Convergence Summary (v{VERSION})",
        "",
        f"**Source:** `{meta['source']}`  ",
        f"**Format:** {meta['source_format']}  ",
        f"**Mode:** {meta['calculation_mode']}  ",
        f"**Ionic Status:** {meta['ionic_status']}  ",
        f"**Ending:** {meta['ending_status']}  ",
        f"**Completed Ionic/Static Evaluations:** {meta['completed_ionic_evaluations']}  ",
        f"**Electronic Cycles Explicitly Confirmed in OUTCAR:** {meta['electronic_cycles_confirmed_by_OUTCAR']}  ",
    ]
    if meta.get("mask_description"):
        lines.append(f"**Force Mask:** {meta['mask_description']}  ")
    lines.append("")
    
    selected = ["EDIFF", "EDIFFG", "NELM", "NSW", "IBRION", "ISIF", "POTIM", "LDIPOL"]
    settings_str = ", ".join(f"{key}={settings[key]}" for key in selected if key in settings)
    lines.append(f"**Settings:** {settings_str}  ")
    if settings.get("EDIFFG", 0) < 0:
        lines.append(f"**Free-atom Force Target:** < {abs(settings['EDIFFG']):.6f} eV/Å  ")
    lines.append("")
    
    lines.append("## Final Ionic Energies and Forces")
    lines.append("*(Not intermediate SCF energies)*")
    lines.append("")
    lines.append("| Eval | SCF* | Final F [eV] | Adjacent dF | Fmax Free | Fmax Raw | Mag [μB] | SCF Status |")
    lines.append("|:-----:|--------:|-------------:|-------------:|----------:|----------:|---------:|------------|")
    
    visible = steps[-last_count:] if last_count else steps
    for step in visible:
        count = str(step["scf_iterations_started"])
        if step["scf_iterations_with_energy"] != step["scf_iterations_started"]:
            count += f"/{step['scf_iterations_with_energy']}"
        unfinished = " [unfinished]" if not step["ionic_complete"] else ""
        lines.append(f"| {step['evaluation']:4d} | {count:>6} | {fmt(step['final_F_eV']):>14} | {fmt(step['delta_F_eV']):>12} | {fmt(step['free_max_force_eV_A']):>11} | {fmt(step['raw_max_force_eV_A']):>10} | {fmt(step['magnetization_muB'], 4):>12} | {step['scf_convergence']}{unfinished} |")
    
    lines.append("")
    lines.append("* SCF = iterations started; started/with-energy shown when the last one is incomplete.")
    lines.append("* Adjacent dF is computed from consecutive final F values; it is not VASP's line-search dE column.")
    
    if steps:
        last = steps[-1]
        lines.append("")
        lines.append(f"## Latest Electronic Cycle (at evaluation {last['evaluation']})")
        lines.append("")
        lines.append(f"- **Iterations started / with an energy:** {last['scf_iterations_started']} / {last['scf_iterations_with_energy']}")
        lines.append(f"- **Last SCF energy:** {fmt(last['last_scf_energy_eV'])} eV; **dE:** {fmt(last['last_scf_dE_eV'], 9)} eV; **d eps:** {fmt(last['last_scf_d_eps_eV'], 9)} eV")
        lines.append("- *Electronic energies can differ from final ionic F (for example, later D3 terms). Compare like with like.*")
    
    lines.append("")
    lines.append("## Review Flags (Heuristics, Not Convergence Criteria)")
    lines.append("")
    if meta["review_flags"]:
        for flag in meta["review_flags"]:
            lines.append(f"- {flag}")
    else:
        lines.append("- No threshold-based escalation flags; this is not a convergence certificate.")
    
    warnings = Counter(item["message"] for item in meta["warnings"])
    lines.append("")
    lines.append("## Warnings / Line-Search Messages")
    lines.append("")
    if warnings:
        for message, count in warnings.items():
            lines.append(f"- {count}x {message}")
    else:
        lines.append("- None recognized in this input snapshot.")
    
    for override in meta.get("mass_overrides", []):
        lines.append(f"- Species {override['species_index']} mass override: INCAR={override['INCAR_mass']:g}, potential metadata in OUTCAR={override['potential_mass']:g}")
    
    lines.append("")
    lines.append("## Errors")
    lines.append("")
    if meta["errors"]:
        for error in meta["errors"]:
            lines.append(f"- ERROR at line {error['line']}: {error['message']}")
    else:
        lines.append("- None")
    
    lines.append("")
    lines.append("## Notes and Limitations")
    lines.append("")
    for note in meta["limitations"]:
        lines.append(f"- {note}")
    for step in steps:
        for note in step["notes"]:
            lines.append(f"- Evaluation {step['evaluation']}: {note}")
    
    lines.append("")
    lines.append("## Official Interpretation References")
    lines.append("")
    lines.append("- [VASP Wiki: OSZICAR](https://vasp.at/wiki/OSZICAR)")
    lines.append("- [VASP Wiki: EDIFF](https://vasp.at/wiki/EDIFF)")
    lines.append("- [VASP Wiki: EDIFFG](https://vasp.at/wiki/EDIFFG)")
    lines.append("- [VASP Wiki: OUTCAR](https://vasp.at/wiki/OUTCAR)")
    
    return "\n".join(lines) + "\n"


def report_text(steps, meta, last_count):
    settings = meta["settings"]
    lines = ["VASP convergence summary (v%s)" % VERSION,
             "Source: %s" % meta["source"], "Format: %s" % meta["source_format"],
             "Mode: %s" % meta["calculation_mode"],
             "Ionic status: %s" % meta["ionic_status"],
             "Ending: %s" % meta["ending_status"],
             "Completed ionic/static evaluations: %d" % meta["completed_ionic_evaluations"],
             "Electronic cycles explicitly confirmed in OUTCAR: %d" % meta["electronic_cycles_confirmed_by_OUTCAR"]]
    if meta.get("mask_description"):
        lines.append("Force mask: " + meta["mask_description"])
    selected = ["EDIFF", "EDIFFG", "NELM", "NSW", "IBRION", "ISIF", "POTIM", "LDIPOL"]
    lines.append("Settings: " + ", ".join("%s=%s" % (key, settings[key]) for key in selected if key in settings))
    if settings.get("EDIFFG", 0) < 0:
        lines.append("Free-atom force target: < %.6f eV/A" % abs(settings["EDIFFG"]))
    lines.extend(["", "FINAL IONIC ENERGIES AND FORCES (not intermediate SCF energies)",
                  "Eval  SCF*  Final F [eV]   adjacent dF   Fmax free   Fmax raw    mag [muB]  SCF status"])
    visible = steps[-last_count:] if last_count else steps
    for step in visible:
        count = str(step["scf_iterations_started"])
        if step["scf_iterations_with_energy"] != step["scf_iterations_started"]:
            count += "/%d" % step["scf_iterations_with_energy"]
        lines.append("%4d %6s %14s %12s %11s %10s %12s  %s%s" % (
            step["evaluation"], count, fmt(step["final_F_eV"]), fmt(step["delta_F_eV"]),
            fmt(step["free_max_force_eV_A"]), fmt(step["raw_max_force_eV_A"]),
            fmt(step["magnetization_muB"], 4), step["scf_convergence"],
            " [unfinished]" if not step["ionic_complete"] else ""))
    lines.append("* SCF = iterations started; started/with-energy shown when the last one is incomplete.")
    lines.append("Adjacent dF is computed from consecutive final F values; it is not VASP's line-search d E column.")
    if steps:
        last = steps[-1]
        lines.extend(["", "LATEST ELECTRONIC CYCLE (at evaluation %d)" % last["evaluation"],
                      "Iterations started / with an energy: %d / %d" % (last["scf_iterations_started"], last["scf_iterations_with_energy"]),
                      "Last SCF energy: %s eV; dE=%s eV; d eps=%s eV" % (fmt(last["last_scf_energy_eV"]), fmt(last["last_scf_dE_eV"], 9), fmt(last["last_scf_d_eps_eV"], 9)),
                      "Electronic energies can differ from final ionic F (for example, later D3 terms). Compare like with like."])
    lines.extend(["", "REVIEW FLAGS (heuristics, not convergence criteria)"])
    lines.extend("- " + flag for flag in meta["review_flags"])
    if not meta["review_flags"]:
        lines.append("- No threshold-based escalation flags; this is not a convergence certificate.")
    warnings = Counter(item["message"] for item in meta["warnings"])
    lines.extend(["", "WARNINGS / LINE-SEARCH MESSAGES"])
    lines.extend("- %dx %s" % (count, message) for message, count in warnings.items())
    for override in meta.get("mass_overrides", []):
        lines.append("- Species %d mass override: INCAR=%g, potential metadata in OUTCAR=%g" % (override["species_index"], override["INCAR_mass"], override["potential_mass"]))
    if not warnings and not meta.get("mass_overrides"):
        lines.append("- None recognized in this input snapshot.")
    lines.extend("- ERROR at line %d: %s" % (x["line"], x["message"]) for x in meta["errors"])
    lines.extend(["", "NOTES AND LIMITATIONS"])
    lines.extend("- " + note for note in meta["limitations"])
    for step in steps:
        lines.extend("- Evaluation %d: %s" % (step["evaluation"], note) for note in step["notes"])
    lines.extend(["", "Official interpretation references:",
                  "https://vasp.at/wiki/OSZICAR", "https://vasp.at/wiki/EDIFF",
                  "https://vasp.at/wiki/EDIFFG", "https://vasp.at/wiki/OUTCAR"])
    return "\n".join(lines) + "\n"


def write_reports(prefix, steps, meta, report, md_report, protected):
    prefix = Path(prefix).expanduser()
    outputs = {"text": Path(str(prefix) + ".txt"), "json": Path(str(prefix) + ".json"),
               "ionic": Path(str(prefix) + "_ionic.csv"), "scf": Path(str(prefix) + "_scf.csv"),
               "md": Path(str(prefix) + ".md")}
    if any(path.resolve() in protected for path in outputs.values()):
        raise ValueError("A report path would overwrite an input; choose another --save prefix")
    for path in outputs.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    outputs["text"].write_text(report, encoding="utf-8")
    outputs["json"].write_text(json.dumps({"metadata": meta, "evaluations": steps}, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    outputs["md"].write_text(md_report, encoding="utf-8")
    columns = ["evaluation", "ionic_complete", "scf_convergence", "scf_iterations_started",
               "scf_iterations_with_energy", "final_F_eV", "delta_F_eV", "E0_eV",
               "last_scf_energy_eV", "last_scf_dE_eV", "last_scf_d_eps_eV", "last_scf_rms_charge",
               "free_max_force_eV_A", "raw_max_force_eV_A", "free_max_force_atom",
               "free_force_below_target", "magnetization_muB", "pressure_kB",
               "force_drift_norm_eV_A", "dispersion_eV"]
    with outputs["ionic"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(steps)
    scf_columns = ["evaluation", "iteration", "energy_eV", "dE_eV", "d_eps_eV", "rms_charge", "magnetization_muB", "loop_seconds", "line"]
    with outputs["scf"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=scf_columns, extrasaction="ignore")
        writer.writeheader()
        for step in steps:
            writer.writerows(step["scf"])
    return outputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("path", nargs="?", default=".", help="Run directory or a log file; defaults to current directory")
    parser.add_argument("--poscar", type=Path, help="Compatible POSCAR with the run's atom order and constraints")
    parser.add_argument("--incar", type=Path, help="INCAR fallback when settings are absent from the log")
    parser.add_argument("--format", choices=["auto", "outcar", "oszicar", "stdout"], default="auto")
    parser.add_argument("--save", type=Path, metavar="PREFIX", help="Write/replace PREFIX.txt, .json, _ionic.csv, _scf.csv; inputs stay read-only")
    parser.add_argument("--md", type=Path, metavar="MD_PREFIX", help="Write/replace MD_PREFIX.md; inputs stay read-only")
    parser.add_argument("--last", type=int, default=20, help="Evaluations to show in text report; 0 means all; CSV/JSON always contain all")
    parser.add_argument("--assume-all-free", action="store_true", help="Explicitly treat all atoms as free when POSCAR is absent (e.g. unconstrained bulk)")
    parser.add_argument("--energy-jump", type=float, default=1.0, help="Review flag for positive adjacent ionic F change in eV/cell (default 1)")
    parser.add_argument("--force-warning", type=float, default=1.0, help="Review flag when free force exceeds this AND triples (eV/A; default 1)")
    parser.add_argument("--drift-warning", type=float, default=1.0, help="Review flag for reported force-drift norm in eV/A (default 1)")
    args = parser.parse_args(argv)
    if args.last < 0 or any(not math.isfinite(x) or x <= 0 for x in [args.energy_jump, args.force_warning, args.drift_warning]):
        parser.error("--last must be nonnegative and review thresholds must be positive finite numbers")
    path = Path(args.path).expanduser()
    if path.is_dir():
        candidates = [path/name for name in ["OUTCAR", "OUTCAR.gz", "OSZICAR", "OSZICAR.gz"]]
        source = next((p for p in candidates if p.is_file()), None)
        if source is None:
            parser.error("No OUTCAR/OSZICAR found; pass a stdout log filename explicitly")
        run_dir = path
    else:
        source, run_dir = path, path.parent
    try:
        text = read_text(source)
        if not text.strip():
            raise ValueError("Log is empty; VASP may not have started")
        source_format = args.format
        if source_format == "auto":
            source_format = "outcar" if re.search(r"Iteration\s+\d+\(\s*\d+\s*\)", text) or "INCAR:" in text else "oszicar"
        poscar_path = args.poscar or (run_dir / "POSCAR" if (run_dir / "POSCAR").is_file() else None)
        incar_path = args.incar or (run_dir / "INCAR" if (run_dir / "INCAR").is_file() else None)
        mask = read_poscar(poscar_path) if poscar_path else None
        if args.assume_all_free and mask:
            raise ValueError("Do not combine --assume-all-free with a POSCAR; its actual flags take priority")
        fallback = read_settings(incar_path)
        meta = {"script_version": VERSION, "source": str(source.resolve()), "source_format": source_format,
                "poscar": str(poscar_path.resolve()) if poscar_path else None,
                "incar_fallback": str(incar_path.resolve()) if incar_path else None,
                "settings": {}, "warnings": [], "errors": [], "limitations": [],
                "normal_footer_seen": False, "structural_accuracy_message": False,
                "stop_request_seen": False, "atom_count": None}
        if source_format == "outcar":
            steps = parse_outcar(text.splitlines(), meta, mask)
            for key, value in fallback.items():
                if key not in meta["settings"]:
                    meta["settings"][key] = value
                    meta["limitations"].append("%s taken from supplied/current INCAR, not confirmed in log" % key)
            if mask and (meta["atom_count"] != mask["atoms"] or (meta.get("ions_per_type") and meta["ions_per_type"] != mask["counts"])):
                raise ValueError("POSCAR counts do not match OUTCAR; use the matching run's POSCAR")
            if args.assume_all_free:
                n = meta["atom_count"]
                if not n:
                    raise ValueError("Cannot assume all free without a known atom count")
                mask = {"atoms": n, "flags": [[True]*3 for _ in range(n)], "free_atoms": n, "fixed_atoms": 0, "partial": False}
                # Re-parse force tables now that the explicitly requested mask is known.
                temp = {**meta, "warnings": [], "errors": []}
                steps = parse_outcar(text.splitlines(), temp, mask)
                meta["mask_description"] = "%d atoms assumed free by explicit option" % n
        else:
            steps = parse_compact(text.splitlines(), meta, fallback)
        if not steps:
            raise ValueError("No supported DFT electronic/ionic records found; a setup-only or other calculation type is not assessable")
        if mask and "mask_description" not in meta:
            meta["mask_description"] = "%d free, %d fixed, %d total from POSCAR%s" % (mask["free_atoms"], mask["fixed_atoms"], mask["atoms"], "; partial directional constraints present" if mask["partial"] else "")
        finish_analysis(steps, meta, args, mask)
        report = report_text(steps, meta, args.last)
        md_report = report_markdown(steps, meta, args.last)
        print(report, end="")
        if args.save:
            protected = {p.resolve() for p in [source, poscar_path, incar_path] if p is not None}
            outputs = write_reports(args.save, steps, meta, report, md_report, protected)
            print("\nReports written:")
            for output in outputs.values():
                print("  " + str(output.resolve()))
        if args.md:
            protected_md = {p.resolve() for p in [source, poscar_path, incar_path] if p is not None}
            md_path = Path(args.md).expanduser()
            if md_path.resolve() in protected_md:
                raise ValueError("MD report path would overwrite an input; choose another --md prefix")
            md_path.parent.mkdir(parents=True, exist_ok=True)
            md_path.write_text(md_report, encoding="utf-8")
            print(f"\nMarkdown report written: {md_path.resolve()}")
        return 0  # Analysis success does not mean the VASP calculation converged.
    except (OSError, ValueError, IndexError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
