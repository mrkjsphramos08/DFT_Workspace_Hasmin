#!/usr/bin/env python3
"""
surface_energy.py

Surface energy of a Quantum ESPRESSO (pw.x) slab calculation:

                              E_slab - n_fu * E_bulk_per_formula_unit
gamma [eV/Ang^2]  =  -----------------------------------------------------
                                  n_surfaces * A_surface

with
    E_slab                   total energy of the slab cell,
    E_bulk_per_formula_unit  bulk total energy / number of formula units in
                             the bulk cell,
    n_fu                     number of bulk formula units in the slab cell,
    n_surfaces               number of equivalent exposed surfaces (default 2: a free-
                             standing slab exposes top and bottom),
    A_surface                in-plane area of the slab cell (Ang^2).

Usage (the default paths point at this project's calculations):
    python surface_energy.py
    python surface_energy.py --slab 04_slab/3_layers/3x2_3_layers.out
    python surface_energy.py --bulk 01_vc_relax/LCO_bulk.out
    python surface_energy.py --n-surfaces 1   # single exposed surface
    python surface_energy.py --area 96.46 --n-fu 18

Everything is read from the .out files themselves:
  * energies  -> the last converged "!    total energy" line of each run;
  * slab cell -> the CELL_PARAMETERS block that belongs to that energy, or, for
                 a plain 'relax' run (which prints no CELL_PARAMETERS), the
                 "crystal axes" block scaled by celldm(1);
  * area      -> |a_i x a_j| of the two in-plane lattice vectors; the vacuum
                 axis is auto-detected, otherwise use --in-plane ab/ac/bc;
  * formula units -> element counts parsed from the atomic positions block,
                 reduced with the bulk composition. Atom labels carrying a
                 numeric suffix (e.g. Co1/Co2, as QE prints when different
                 sublattices are given different starting_magnetization
                 under nspin=2) are folded back into the base element so
                 they are not silently dropped from the stoichiometry;
  * surfaces  -> number of equivalent exposed surfaces (default 2 for a free-
                 standing slab: top and bottom), override with --n-surfaces.

The script also verifies that slab and bulk share the same settings (pseudo-
potentials, cutoffs, functional, U, dispersion, smearing) and flags relax /
vc-relax runs that are still in progress, because in that case the reference
energies are only provisional.
"""

import argparse
import re
import sys
from functools import reduce
from math import gcd
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Constants and default inputs
# ---------------------------------------------------------------------------
# Floating point number as printed by QE, Fortran 'D' exponents allowed.
NUM = r"([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)"

RY_TO_EV = 13.605698065894      # 1 Ry in eV
EV_A2_TO_J_M2 = 16.0218         # 1 eV/Ang^2 in J/m^2
BOHR_TO_ANG = 0.529177249       # 1 Bohr in Angstrom

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SLAB = SCRIPT_DIR / "04_slab" / "4_layers" / "3x2_4_layers.out"
DEFAULT_BULK = SCRIPT_DIR / "01_vc_relax" / "LCO_bulk.out"

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------
def to_float(token):
    """Convert a number printed by QE (possibly with a Fortran 'D' exponent)."""
    return float(str(token).replace("D", "E").replace("d", "e"))


def read_lines(path):
    """All lines of a text file, or exit with a clear message when it is missing."""
    path = Path(path)
    if not path.is_file():
        sys.exit(f"ERROR: file not found: {path}")
    return path.read_text(errors="ignore").splitlines()


def parse_alat_ang(lines):
    """alat = celldm(1) in Angstrom (None when it is not printed)."""
    for line in lines:
        m = re.search(r"celldm\(1\)\s*=\s*" + NUM, line)
        if m:
            return to_float(m.group(1)) * BOHR_TO_ANG
    for line in lines:
        m = re.search(r"lattice parameter \(alat\)\s*=\s*" + NUM, line)
        if m:
            return to_float(m.group(1)) * BOHR_TO_ANG
    return None


def parse_cell_blocks(lines, alat_ang):
    """Every CELL_PARAMETERS block as (line index, 3x3 cell in Angstrom)."""
    blocks = []
    for i, line in enumerate(lines):
        if not line.strip().upper().startswith("CELL_PARAMETERS"):
            continue
        header = line.split("CELL_PARAMETERS", 1)[1].strip().strip("(){}[]").lower()
        if "angstrom" in header:
            scale = 1.0
        elif "bohr" in header:
            scale = BOHR_TO_ANG
        else:                                   # 'alat': QE prints alat=<Bohr>
            m = re.search(NUM, header)
            alat = to_float(m.group(1)) * BOHR_TO_ANG if m else alat_ang
            if alat is None:
                continue
            scale = alat

        vecs = []
        for j in (1, 2, 3):
            parts = lines[i + j].split() if i + j < len(lines) else []
            if len(parts) < 3:
                vecs = None
                break
            try:
                vecs.append([to_float(x) for x in parts[:3]])
            except ValueError:
                vecs = None
                break
        if vecs is not None:
            blocks.append((i, np.array(vecs) * scale))
    return blocks


def parse_crystal_axes(lines, alat_ang):
    """The 3x3 cell (Angstrom) from the 'crystal axes' block (units of alat)."""
    if alat_ang is None:
        return None
    rows = {}
    for line in lines:
        m = re.match(r"\s*a\((\d)\)\s*=\s*\(\s*" + NUM + r"\s+" + NUM + r"\s+" + NUM + r"\s*\)", line)
        if m:
            rows[int(m.group(1))] = [to_float(m.group(k)) for k in (2, 3, 4)]
    if len(rows) == 3:
        return np.array([rows[1], rows[2], rows[3]]) * alat_ang
    return None


def parse_converged_energies(lines):
    """[(line index, energy Ry)] for every converged '!    total energy' line."""
    pattern = re.compile(r"^\s*!\s*total energy\s*=\s*" + NUM + r"\s*Ry")
    found = []
    for i, line in enumerate(lines):
        m = pattern.match(line)
        if m:
            found.append((i, to_float(m.group(1))))
    return found


def parse_species_counts(lines):
    """{base_element: number of atoms} from the FIRST 'positions (alat units)'
    block in the file.

    QE reprints this atomic-summary table every time the cell changes -
    e.g. at each accepted step of a vc-relax run - so scanning the whole file
    and counting every matching line multiplies the true atom count by
    however many times the block happens to be reprinted. The composition
    itself does not change between those reprints, so only the first block
    is parsed.

    QE atom labels may carry a numeric suffix (e.g. Co1 / Co2) when different
    sublattices are assigned different starting_magnetization, which is the
    normal way to set up an antiferromagnetic or per-site Hubbard U treatment
    under nspin=2. Those suffixed labels are folded back into their base
    element (Co1, Co2 -> Co) so the stoichiometry doesn't silently lose atoms.
    """
    header_re = re.compile(r"positions\s*\(alat units\)", re.IGNORECASE)
    atom_re = re.compile(r"^\s*\d+\s+([A-Za-z]{1,3}\d*)\s+tau\(")
    counts = {}
    i, n = 0, len(lines)
    while i < n:
        if header_re.search(lines[i]):
            i += 1
            while i < n:
                m = atom_re.match(lines[i])
                if not m:
                    break
                element = re.sub(r"\d+$", "", m.group(1))
                counts[element] = counts.get(element, 0) + 1
                i += 1
            break  # only the first block is used - composition is static
        i += 1
    return counts


def parse_n_atoms(lines):
    """Value of 'number of atoms/cell' (None when it is not printed)."""
    for line in lines:
        m = re.search(r"number of atoms/cell\s*=\s*(\d+)", line)
        if m:
            return int(m.group(1))
    return None


def parse_status(lines):
    """Run and geometry-optimisation status flags of a pw.x output."""
    return {
        "finished": any(line.strip().startswith("JOB DONE") for line in lines),
        "is_relax": any("BFGS Geometry Optimization" in line for line in lines),
        "relax_converged": any("bfgs converged" in line for line in lines),
    }


def parse_settings(lines):
    """The pw.x settings that have to agree between the slab and bulk runs."""
    text = "\n".join(lines)

    def number(pattern):
        m = re.search(pattern, text)
        return to_float(m.group(1)) if m else None

    settings = {
        "ecutwfc": number(r"kinetic-energy cutoff\s*=\s*" + NUM + r"\s*Ry"),
        "ecutrho": number(r"charge density cutoff\s*=\s*" + NUM + r"\s*Ry"),
        "dispersion": "DFT-D3" if "DFT-D3 Dispersion Correction" in text else "none",
        "spin_polarized": ("Starting magnetic structure" in text) or ("SPIN UP" in text),
    }

    m = re.search(r"Exchange-correlation=\s*(.+)", text)
    settings["xc"] = " ".join(m.group(1).split()) if m else None

    m = re.search(r"(\w+) smearing, width \(Ry\)=\s*" + NUM, text)
    settings["smearing"] = (m.group(1), to_float(m.group(2))) if m else None

    settings["hubbard"] = {m.group(1): to_float(m.group(2))
                           for m in re.finditer(r"U\(([^)]+)\)\s*=\s*" + NUM, text)}

    settings["pseudos"] = {
        m.group(1): Path(m.group(2)).name
        for m in re.finditer(
            r"PseudoPot\.\s*#\s*\d+\s+for\s+(\w+)\s+read from file:\s*\n\s*(\S+)", text)
    }

    return settings


def load_run(path):
    """Parse everything the surface-energy calculation needs from one .out file."""
    lines = read_lines(path)
    alat_ang = parse_alat_ang(lines)
    return {
        "path": Path(path),
        "lines": lines,
        "alat_ang": alat_ang,
        "cell_blocks": parse_cell_blocks(lines, alat_ang),
        "crystal_axes": parse_crystal_axes(lines, alat_ang),
        "energies": parse_converged_energies(lines),
        "n_atoms": parse_n_atoms(lines),
        "counts": parse_species_counts(lines),
        "status": parse_status(lines),
        "settings": parse_settings(lines),
    }

# ---------------------------------------------------------------------------
# Structure and energy helpers
# ---------------------------------------------------------------------------
def reference_energy(run):
    """Last converged SCF energy (Ry) and the cell used for that energy.

    For vc-relax runs every converged energy is paired with the
    CELL_PARAMETERS block that precedes it, so a trailing cell block without a
    following converged energy (run still in progress) is ignored.
    """
    if not run["energies"]:
        sys.exit(f"ERROR: no converged '!    total energy' line found in {run['path']}")
    line_no, energy = run["energies"][-1]
    cell = None
    for start, block in run["cell_blocks"]:
        if start < line_no:
            cell = block
    if cell is None:
        cell = run["crystal_axes"]
    return energy, cell


def detect_in_plane_axes(cell):
    """Indices of the two lattice vectors that span the surface plane.

    The vacuum axis is the vector whose perpendicular distance to the plane of
    the other two vectors is the largest; those other two are the in-plane axes.
    """
    vacuum_axis, largest = 2, -1.0
    for k in range(3):
        others = [cell[m] for m in range(3) if m != k]
        normal = np.cross(others[0], others[1])
        norm = float(np.linalg.norm(normal))
        if norm < 1e-8:
            continue
        distance = abs(float(np.dot(cell[k], normal))) / norm
        if distance > largest:
            largest, vacuum_axis = distance, k
    in_plane = tuple(m for m in range(3) if m != vacuum_axis)
    return in_plane, vacuum_axis


def surface_area(cell, in_plane=None):
    """Area of the in-plane cell in Angstrom^2."""
    if in_plane is None:
        in_plane, _ = detect_in_plane_axes(cell)
    i, j = in_plane
    return float(np.linalg.norm(np.cross(cell[i], cell[j])))


def reduced_formula(counts):
    """Reduce an element-count dictionary to its primitive formula unit."""
    values = [c for c in counts.values() if c]
    if not values:
        return {}, 0
    divisor = reduce(gcd, values)
    return {element: c // divisor for element, c in counts.items()}, divisor


def formula_units(counts, fu_counts):
    """Number n such that counts == n * fu_counts, or None when impossible."""
    per_fu = sum(fu_counts.values())
    if per_fu == 0 or set(counts) != set(fu_counts):
        return None
    n = sum(counts.values()) / per_fu
    if abs(n - round(n)) > 1e-6:
        return None
    n = int(round(n))
    if any(counts[element] != n * c for element, c in fu_counts.items()):
        return None
    return n


def format_formula(counts):
    """'Li18 Co18 O36' style string for an element-count dictionary."""
    return " ".join(f"{element}{counts[element]}" for element in sorted(counts))

# ---------------------------------------------------------------------------
# Reporting and command line interface
# ---------------------------------------------------------------------------
def compare_settings(slab_settings, bulk_settings):
    """Rows (label, slab value, bulk value, identical) for the shared settings."""
    rows = []
    keys = (
        ("ecutwfc", "ecutwfc [Ry]"),
        ("ecutrho", "ecutrho [Ry]"),
        ("xc", "XC functional"),
        ("dispersion", "dispersion"),
        ("smearing", "smearing"),
        ("spin_polarized", "spin polarised"),
    )
    for key, label in keys:
        slab_value = slab_settings.get(key)
        bulk_value = bulk_settings.get(key)
        rows.append((label, slab_value, bulk_value, slab_value == bulk_value))

    for element in sorted(set(slab_settings["pseudos"]) | set(bulk_settings["pseudos"])):
        slab_value = slab_settings["pseudos"].get(element)
        bulk_value = bulk_settings["pseudos"].get(element)
        rows.append((f"pseudo {element}", slab_value, bulk_value, slab_value == bulk_value))

    for element in sorted(set(slab_settings["hubbard"]) | set(bulk_settings["hubbard"])):
        slab_value = slab_settings["hubbard"].get(element)
        bulk_value = bulk_settings["hubbard"].get(element)
        rows.append((f"U({element})", slab_value, bulk_value, slab_value == bulk_value))

    return rows


def format_value(value):
    """Compact representation of a setting value."""
    if isinstance(value, float):
        return f"{value:g}"
    if isinstance(value, tuple):
        return " ".join(format_value(item) for item in value)
    return str(value)


def run_status_text(run):
    """Human readable status of one run, flagging output that is incomplete."""
    status = run["status"]
    if not status["finished"]:
        return "unfinished - no 'JOB DONE.', energies are provisional"
    if status["is_relax"] and not status["relax_converged"]:
        return "unfinished - geometry optimisation not converged"
    if status["is_relax"]:
        return "ok - geometry optimisation converged"
    return "ok"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Surface energy from Quantum ESPRESSO slab and bulk outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--slab", default=str(DEFAULT_SLAB), help="slab .out file")
    parser.add_argument("--bulk", default=str(DEFAULT_BULK), help="bulk .out file")
    parser.add_argument("--n-fu", type=float, default=None,
                        help="bulk formula units in the slab cell (default: auto)")
    parser.add_argument("--area", type=float, default=None,
                        help="surface area in Ang^2 (default: from the slab cell)")
    parser.add_argument("--in-plane", default=None, choices=("ab", "ac", "bc"),
                        help="in-plane lattice vectors (default: auto)")
    parser.add_argument("--n-surfaces", type=int, default=2, metavar="N",
                        help="number of equivalent exposed surfaces; 2 for a "
                             "free-standing slab (top and bottom), 1 for a "
                             "single exposed surface")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    slab = load_run(args.slab)
    bulk = load_run(args.bulk)

    # --- formula units in the slab and in the bulk cell ---------------------
    fu_counts, _ = reduced_formula(bulk["counts"])
    if not fu_counts:
        sys.exit(f"ERROR: no atomic species parsed from {bulk['path']}.")

    n_fu_bulk = formula_units(bulk["counts"], fu_counts)
    if n_fu_bulk is None:
        sys.exit(f"ERROR: cannot determine the number of formula units in "
                 f"{bulk['path']} (composition {format_formula(bulk['counts'])}).")

    n_fu_slab_auto = formula_units(slab["counts"], fu_counts)
    if args.n_fu is not None:
        n_fu_slab = args.n_fu
    else:
        if n_fu_slab_auto is None:
            sys.exit(f"ERROR: slab composition {format_formula(slab['counts'])} is not an "
                     f"integer multiple of the bulk formula unit {format_formula(fu_counts)}; "
                     f"pass --n-fu to define the number of formula units.")
        n_fu_slab = n_fu_slab_auto

    # --- reference energies and slab cell -----------------------------------
    e_slab_ry, slab_cell = reference_energy(slab)
    e_bulk_ry, _ = reference_energy(bulk)
    e_bulk_fu_ry = e_bulk_ry / n_fu_bulk

    # --- surface area -------------------------------------------------------
    in_plane = tuple("abc".index(char) for char in args.in_plane) if args.in_plane else None
    if args.area is not None:
        area = args.area
        area_note = "given on the command line"
    elif slab_cell is not None:
        if in_plane is None:
            in_plane, _ = detect_in_plane_axes(slab_cell)
        area = surface_area(slab_cell, in_plane)
        area_note = f"|a{in_plane[0] + 1} x a{in_plane[1] + 1}| of the slab cell"
    else:
        sys.exit(f"ERROR: no cell information found in {slab['path']}; pass --area.")

    # --- number of equivalent exposed surfaces ------------------------------
    # A free-standing slab exposes two surfaces (top and bottom), so the
    # default is n_surfaces = 2. Use --n-surfaces 1 for a single exposed
    # surface (e.g. one side is a substrate, or only one surface is of
    # interest).
    n_surfaces = args.n_surfaces
    surfaces_note = "two (top and bottom)" if n_surfaces == 2 else str(n_surfaces)

    # --- surface energy -----------------------------------------------------
    numerator_ry = e_slab_ry - n_fu_slab * e_bulk_fu_ry
    gamma_ev_a2 = numerator_ry * RY_TO_EV / (n_surfaces * area)
    gamma_j_m2 = gamma_ev_a2 * EV_A2_TO_J_M2

    def row(label, value):
        print(f"  {label:<28}: {value}")

    print("=" * 78)
    print("QE slab surface energy")
    print("=" * 78)
    print(f"  SLAB  {slab['path']}")
    row("run status", run_status_text(slab))
    row("atoms / cell", slab["n_atoms"])
    row("composition", f"{format_formula(slab['counts'])}"
                       f"  =  {n_fu_slab:g} x {format_formula(fu_counts)}")
    row("E_slab (F)", f"{e_slab_ry:.8f} Ry  =  {e_slab_ry * RY_TO_EV:.4f} eV")
    row("surface area", f"{area:.4f} Ang^2   ({area_note})")
    row("exposed surfaces", f"{n_surfaces}   ({surfaces_note})")
    print()
    print(f"  BULK  {bulk['path']}")
    row("run status", run_status_text(bulk))
    row("atoms / cell", bulk["n_atoms"])
    row("composition", f"{format_formula(bulk['counts'])}"
                       f"  =  {n_fu_bulk} x {format_formula(fu_counts)}")
    row("E_bulk (F)", f"{e_bulk_ry:.8f} Ry  =  {e_bulk_ry * RY_TO_EV:.4f} eV")
    row("E_bulk / formula unit", f"{e_bulk_fu_ry:.8f} Ry  =  {e_bulk_fu_ry * RY_TO_EV:.4f} eV")
    print()
    print("  SETTINGS CONSISTENCY (slab vs bulk)")

    settings_rows = compare_settings(slab["settings"], bulk["settings"])
    for label, slab_value, bulk_value, identical in settings_rows:
        mark = "ok" if identical else "DIFFERS"
        row(label, f"{format_value(slab_value)}  vs  {format_value(bulk_value)}   [{mark}]")

    print()
    print("  RESULT")
    row("E_slab - n_fu * E_bulk/f.u.", f"{numerator_ry:.8f} Ry"
                                      f"  =  {numerator_ry * RY_TO_EV:.4f} eV")
    row("n_surfaces * area", f"{n_surfaces * area:.4f} Ang^2")
    row("surface energy", f"{gamma_ev_a2:.6f} eV/Ang^2"
                          f"  =  {gamma_j_m2:.4f} J/m^2")
    print("=" * 78)

    # --- warnings -----------------------------------------------------------
    warnings = []
    for run in (slab, bulk):
        if not run["status"]["finished"]:
            warnings.append(f"{run['path']}: the run is still in progress (no 'JOB DONE.'); "
                            f"the last converged energy is used as it is.")
        elif run["status"]["is_relax"] and not run["status"]["relax_converged"]:
            warnings.append(f"{run['path']}: the geometry optimisation did not converge.")
    for run in (slab, bulk):
        energies = [energy for _, energy in run["energies"]]
        if len(energies) >= 2:
            jump = abs(energies[-1] - energies[-2])
            if jump > 0.01:
                warnings.append(f"{run['path']}: the last two converged SCF energies differ by "
                                f"{jump:.4f} Ry ({energies[-2]:.6f} -> {energies[-1]:.6f} Ry); the "
                                f"relaxation is jumping between electronic solutions, so this "
                                f"reference energy is not stable - rerun to convergence first.")
    for label, slab_value, bulk_value, identical in settings_rows:
        if not identical:
            warnings.append(f"slab and bulk differ for {label}: "
                            f"{format_value(slab_value)} vs {format_value(bulk_value)}")
    if args.n_fu is not None and n_fu_slab_auto is None:
        warnings.append(f"--n-fu {n_fu_slab:g} was used although the slab composition "
                        f"{format_formula(slab['counts'])} is not a multiple of "
                        f"{format_formula(fu_counts)}.")
    if warnings:
        print("WARNINGS")
        for warning in warnings:
            print(f"  ! {warning}")
        print("=" * 78)

    return {"surface_energy_ev_a2": gamma_ev_a2,
            "surface_energy_j_m2": gamma_j_m2,
            "area": area,
            "n_fu_slab": n_fu_slab,
            "n_surfaces": n_surfaces}


if __name__ == "__main__":
    main()