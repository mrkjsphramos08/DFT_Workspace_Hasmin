#!/usr/bin/env python3
"""Calculate a stoichiometric VASP slab surface energy from OUTCAR files.

The default result is the average energy of the two faces of a slab:
    gamma_avg = (E_slab - N * E_bulk) / (2 * A)

For a slab with only one face relaxed, supply --unrelaxed-slab-outcar to
obtain the common approximation gamma_relaxed = 2*gamma_avg-gamma_unrelaxed.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple


VERSION = "1.0.0"
EV_PER_A2_TO_J_PER_M2 = 16.02176634
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


class ParseError(ValueError):
    """Raised when an OUTCAR does not contain the required final data."""


def read_text(path: Path) -> str:
    try:
        if path.suffix.lower() == ".gz":
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                return handle.read()
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ParseError("Could not read {}: {}".format(path, exc)) from exc


def parse_final_energy(text: str, kind: str, label: str) -> float:
    patterns = {
        "toten": r"free\s+energy\s+TOTEN\s*=\s*(%s)" % FLOAT,
        "without-entropy": r"energy\s+without\s+entropy\s*=\s*(%s)" % FLOAT,
        "e0": r"energy\s*\(\s*sigma\s*->\s*0\s*\)\s*=\s*(%s)" % FLOAT,
    }
    matches = re.findall(patterns[kind], text, flags=re.IGNORECASE)
    if not matches:
        raise ParseError("{} has no final {} energy in the expected OUTCAR format.".format(label, kind))
    return float(matches[-1])


def parse_ion_counts(text: str, label: str) -> Tuple[int, ...]:
    matches = re.findall(r"ions\s+per\s+type\s*=\s*([^\r\n]+)", text, flags=re.IGNORECASE)
    if not matches:
        raise ParseError("{} has no 'ions per type' line.".format(label))
    counts = tuple(int(value) for value in re.findall(r"\d+", matches[-1]))
    if not counts or any(value <= 0 for value in counts):
        raise ParseError("{} has invalid ion counts: {!r}".format(label, counts))
    return counts


def parse_area(text: str, label: str) -> float:
    markers = list(re.finditer(r"direct\s+lattice\s+vectors", text, flags=re.IGNORECASE))
    if not markers:
        raise ParseError("{} has no 'direct lattice vectors' section.".format(label))
    remainder = text[markers[-1].end() :].splitlines()
    vectors: List[Tuple[float, float, float]] = []
    for line in remainder:
        values = re.findall(FLOAT, line)
        if len(values) >= 3:
            vectors.append((float(values[0]), float(values[1]), float(values[2])))
        if len(vectors) == 3:
            break
    if len(vectors) != 3:
        raise ParseError("{} ends before its lattice vectors are complete.".format(label))
    a, b = vectors[0], vectors[1]
    cross = (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )
    area = math.sqrt(sum(value * value for value in cross))
    if area <= 0:
        raise ParseError("{} has zero in-plane area.".format(label))
    return area


def formula_units(counts: Sequence[int]) -> Tuple[int, Tuple[int, ...]]:
    divisor = 0
    for count in counts:
        divisor = math.gcd(divisor, count)
    return divisor, tuple(count // divisor for count in counts)


def gamma(energy_slab: float, energy_bulk: float, multiplier: int, area: float, surfaces: int) -> Tuple[float, float, float]:
    excess = energy_slab - multiplier * energy_bulk
    value_ev_a2 = excess / (surfaces * area)
    return excess, value_ev_a2, value_ev_a2 * EV_PER_A2_TO_J_PER_M2


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate a stoichiometric VASP slab surface energy from final OUTCAR energies."
    )
    parser.add_argument("slab_outcar", type=Path, help="Relaxed slab OUTCAR (or OUTCAR.gz)")
    parser.add_argument("bulk_outcar", type=Path, help="Matching bulk OUTCAR (or OUTCAR.gz)")
    parser.add_argument(
        "--energy",
        choices=("e0", "without-entropy", "toten"),
        default="e0",
        help="Final energy field to use consistently in both files (default: e0, energy(sigma->0)).",
    )
    parser.add_argument(
        "--area",
        type=float,
        metavar="ANGSTROM2",
        help="Override the slab |a x b| area in Å^2 (normally read from the slab OUTCAR).",
    )
    parser.add_argument(
        "--surfaces",
        type=int,
        default=2,
        help="Number of equivalent faces represented by the slab energy (default: 2).",
    )
    parser.add_argument(
        "--unrelaxed-slab-outcar",
        type=Path,
        help="Unrelaxed symmetric slab OUTCAR; reports a top-face relaxed estimate as well.",
    )
    parser.add_argument("--version", action="version", version="surface_energy.py " + VERSION)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    if args.surfaces <= 0:
        raise ParseError("--surfaces must be a positive integer.")
    if args.area is not None and args.area <= 0:
        raise ParseError("--area must be positive.")
    if args.unrelaxed_slab_outcar and args.surfaces != 2:
        raise ParseError("The one-face-relaxed estimate requires --surfaces 2.")

    slab_text = read_text(args.slab_outcar)
    bulk_text = read_text(args.bulk_outcar)
    slab_energy = parse_final_energy(slab_text, args.energy, str(args.slab_outcar))
    bulk_energy = parse_final_energy(bulk_text, args.energy, str(args.bulk_outcar))
    slab_counts = parse_ion_counts(slab_text, str(args.slab_outcar))
    bulk_counts = parse_ion_counts(bulk_text, str(args.bulk_outcar))
    slab_fu, slab_ratio = formula_units(slab_counts)
    bulk_fu, bulk_ratio = formula_units(bulk_counts)
    if len(slab_counts) != len(bulk_counts) or slab_ratio != bulk_ratio:
        raise ParseError(
            "Slab and bulk are not stoichiometrically compatible: slab counts {} (ratio {}), "
            "bulk counts {} (ratio {}). This simple formula is not valid; use chemical potentials."
            .format(slab_counts, slab_ratio, bulk_counts, bulk_ratio)
        )
    if slab_fu % bulk_fu:
        raise ParseError("Slab formula-unit count is not an integer multiple of the bulk cell.")
    multiplier = slab_fu // bulk_fu
    area = args.area if args.area is not None else parse_area(slab_text, str(args.slab_outcar))
    excess, gamma_ev_a2, gamma_j_m2 = gamma(slab_energy, bulk_energy, multiplier, area, args.surfaces)

    print("VASP surface-energy summary (v{})".format(VERSION))
    print("Energy field: {} (used consistently for slab and bulk)".format(args.energy))
    print("Slab: {}".format(args.slab_outcar))
    print("Bulk: {}".format(args.bulk_outcar))
    print("Slab ions per type: {}; inferred formula units: {}".format(" ".join(map(str, slab_counts)), slab_fu))
    print("Bulk ions per type: {}; inferred formula units: {}".format(" ".join(map(str, bulk_counts)), bulk_fu))
    print("Bulk-cell multiplier N: {}".format(multiplier))
    print("Surface area A: {:.8f} A^2".format(area))
    print("E_slab: {:.10f} eV".format(slab_energy))
    print("E_bulk: {:.10f} eV per bulk cell".format(bulk_energy))
    print("Excess energy E_slab - N*E_bulk: {:.10f} eV".format(excess))
    print("Surface energy ({} faces): {:.8f} eV/A^2 = {:.6f} J/m^2".format(args.surfaces, gamma_ev_a2, gamma_j_m2))

    if args.unrelaxed_slab_outcar:
        unrelaxed_text = read_text(args.unrelaxed_slab_outcar)
        unrelaxed_energy = parse_final_energy(unrelaxed_text, args.energy, str(args.unrelaxed_slab_outcar))
        unrelaxed_counts = parse_ion_counts(unrelaxed_text, str(args.unrelaxed_slab_outcar))
        unrelaxed_fu, unrelaxed_ratio = formula_units(unrelaxed_counts)
        if unrelaxed_ratio != slab_ratio or unrelaxed_fu != slab_fu:
            raise ParseError("The unrelaxed slab composition does not match the relaxed slab.")
        unrelaxed_area = parse_area(unrelaxed_text, str(args.unrelaxed_slab_outcar))
        if not math.isclose(unrelaxed_area, area, rel_tol=1e-7, abs_tol=1e-7):
            raise ParseError("The unrelaxed and relaxed slab in-plane areas differ; do not apply the one-face estimate.")
        _, unrelaxed_ev_a2, unrelaxed_j_m2 = gamma(unrelaxed_energy, bulk_energy, multiplier, area, 2)
        relaxed_ev_a2 = 2 * gamma_ev_a2 - unrelaxed_ev_a2
        relaxed_j_m2 = relaxed_ev_a2 * EV_PER_A2_TO_J_PER_M2
        print("Unrelaxed slab: {}".format(args.unrelaxed_slab_outcar))
        print("gamma_unrelaxed: {:.8f} eV/A^2 = {:.6f} J/m^2".format(unrelaxed_ev_a2, unrelaxed_j_m2))
        print("Estimated relaxed top-face gamma: {:.8f} eV/A^2 = {:.6f} J/m^2".format(relaxed_ev_a2, relaxed_j_m2))
        print("This estimate assumes two initially equivalent faces and only one was relaxed.")

    print("\nInterpretation: numerical result only. Verify force, k-point, vacuum, and slab-model convergence separately.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ParseError as exc:
        print("surface_energy.py: error: {}".format(exc), file=sys.stderr)
        raise SystemExit(2)
