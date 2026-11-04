#!/usr/bin/env python3
"""Build and audit a VASP-ready LiCoO2 (104) slab from relaxed bulk CONTCAR.

Dependencies: numpy and pymatgen. No licensed POTCAR data is read or written.

DEFAULT INTERPRETATION: --slab-size 3 --thickness-mode planes requests
three nominal (104) interplanar spacings in pymatgen, NOT three complete
LiCoO2 chemical layers. Inspect the atomic-plane inventory, termination,
stoichiometry, and surface chemistry before using the generated POSCAR.

Example (on Hasmin, using the existing relaxed bulk):
    python lco_104_slab_vasp.py \
      --bulk /home/mark/LCO/01_bulk/01_bulk_vasp/CONTCAR \
      --output-dir ./lco_104_candidate \
      --slab-size 3 --thickness-mode planes \
      --vacuum 18 --vacuum-mode top --supercell 3 2 1 --termination-index 0 \
      --freeze-bottom-planes 2

To reproduce the original script's 6.5 A minimum slab size instead:
    python lco_104_slab_vasp.py --thickness-mode angstrom --slab-size 6.5

Outputs: POSCAR (with Selective dynamics), slab_104.cif,
         slab_report.json, slab_report.txt, and optionally QE_blocks.in.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
import json
import math
from pathlib import Path
import sys

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar

SPECIES_ORDER = ("Li", "Co", "O")
DEFAULT_BULK = "/home/mark/LCO/01_bulk/01_bulk_vasp/CONTCAR"


def gauss_reduce_2d(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unimodular Gauss reduction of two surface-cell vectors."""
    a = np.asarray(a, dtype=float).copy()
    b = np.asarray(b, dtype=float).copy()
    for _ in range(1000):
        if np.dot(b, b) < np.dot(a, a) - 1e-10:
            a, b = b, a
        if np.dot(a, a) < 1e-12:
            raise ValueError("Degenerate surface lattice vector")
        mu = int(np.rint(np.dot(a, b) / np.dot(a, a)))
        if mu == 0:
            break
        b = b - mu * a
    else:
        raise RuntimeError("In-plane Gauss reduction failed to terminate")
    if np.dot(a, b) < 0:
        a = -a
    return a, b


def compact_in_plane(structure: Structure) -> Structure:
    """Shorten in-plane basis without modifying Cartesian atomic geometry."""
    old = np.array(structure.lattice.matrix, dtype=float)
    a, b = gauss_reduce_2d(old[0], old[1])
    new = np.array([a, b, old[2]])
    transform = np.rint(new @ np.linalg.inv(old)).astype(int)

    if not np.allclose(transform @ old, new, atol=1e-7):
        raise RuntimeError("Surface-cell basis change is not integral")
    det = int(round(np.linalg.det(transform)))
    if abs(det) != 1:
        raise RuntimeError(f"Surface-cell basis change has determinant {det}, not +/-1")
    if det < 0:
        new[0] *= -1
        transform = np.rint(new @ np.linalg.inv(old)).astype(int)
    if int(round(np.linalg.det(transform))) != 1:
        raise RuntimeError("Could not obtain a right-handed reduced surface basis")

    new_frac = structure.frac_coords @ np.linalg.inv(transform)
    return Structure(Lattice(new), structure.species, new_frac)


def orient_surface_xy(structure: Structure) -> Structure:
    """Rigidly rotate surface normal onto +z and first cell vector onto +x.

    Requires c to be perpendicular to a,b; never flattens a genuinely
    sheared lattice by deleting vector components (which would strain bonds).
    """
    mat = np.asarray(structure.lattice.matrix, dtype=float)
    ez = mat[2] / np.linalg.norm(mat[2])
    ax = mat[0] - np.dot(mat[0], ez) * ez
    if np.linalg.norm(ax) < 1e-8:
        raise RuntimeError("In-plane a vector is degenerate")
    ex = ax / np.linalg.norm(ax)
    ey = np.cross(ez, ex)
    rotation = np.column_stack([ex, ey, ez])
    new = mat @ rotation
    if not np.allclose(new[:2, 2], 0.0, atol=1e-6):
        raise RuntimeError("c is not perpendicular to surface plane")
    if not np.allclose(new[2, :2], 0.0, atol=1e-6):
        raise RuntimeError("Vacuum direction is not along z")
    if np.linalg.det(new) <= 0:
        raise RuntimeError("Rotated slab cell is left-handed")
    new[np.abs(new) < 1e-12] = 0.0
    return Structure(Lattice(new), structure.species, structure.frac_coords)


def surface_normal_heights(structure: Structure) -> np.ndarray:
    """Coordinates normal to the a-b plane, valid for an initially sheared slab."""
    vectors = np.asarray(structure.lattice.matrix, dtype=float)
    normal = np.cross(vectors[0], vectors[1])
    normal /= np.linalg.norm(normal)
    if np.dot(normal, vectors[2]) < 0:
        normal = -normal
    return np.asarray(structure.cart_coords @ normal, dtype=float)


def group_atomic_planes(z: np.ndarray, tolerance: float) -> list[list[int]]:
    """Cluster atoms by geometric z; groups are NOT chemical LiCoO2 layers."""
    ordered = np.argsort(z)
    groups: list[list[int]] = []
    for idx in ordered:
        i = int(idx)
        if not groups or abs(float(z[i]) - float(np.mean(z[groups[-1]]))) > tolerance:
            groups.append([i])
        else:
            groups[-1].append(i)
    return groups


def species_counts(structure: Structure, indices: list[int] | None = None) -> dict[str, int]:
    if indices is None:
        indices = list(range(len(structure)))
    counts = Counter(structure[i].specie.symbol for i in indices)
    return {s: int(counts[s]) for s in SPECIES_ORDER if counts[s]}


def termination_inventory(slabs: list[Structure], tolerance: float) -> list[dict]:
    info = []
    for i, slab in enumerate(slabs):
        h = surface_normal_heights(slab)
        lower = np.where(h <= h.min() + tolerance)[0].tolist()
        upper = np.where(h >= h.max() - tolerance)[0].tolist()
        info.append({
            "termination_index": i,
            "pymatgen_shift": float(getattr(slab, "shift", 0.0)),
            "atoms_in_unexpanded_cell": len(slab),
            "composition": species_counts(slab),
            "bottom_plane_species": species_counts(slab, lower),
            "top_plane_species": species_counts(slab, upper),
            "slab_atomic_span_angstrom": float(np.ptp(h)),
        })
    return info


def translate_li_xy_origin(structure: Structure) -> tuple[Structure, int]:
    """Put the lowest-z Li atom at x=y=0 modulo periodic lattice translation."""
    li_indices = [i for i, site in enumerate(structure) if site.specie.symbol == "Li"]
    if not li_indices:
        raise RuntimeError("No Li atoms in slab; cannot align Li origin")
    z = np.asarray(structure.cart_coords[:, 2], dtype=float)
    lowest = min(z[i] for i in li_indices)
    near_bottom = [i for i in li_indices if abs(z[i] - lowest) < 1e-3]
    frac = np.array(structure.frac_coords, dtype=float)
    chosen = min(near_bottom, key=lambda i: float(np.linalg.norm(frac[i, :2] % 1.0)))
    frac[:, :2] = (frac[:, :2] - frac[chosen, :2]) % 1.0
    return Structure(structure.lattice, structure.species, frac), chosen


def recenter_with_vacuum(
    structure: Structure,
    vacuum_angstrom: float,
    mode: str = "top",
) -> Structure:
    """Position atoms along z; set geometric vacuum gap exactly to the target.

    If mode == 'top', puts bottom atoms at z=0 and all vacuum above the slab.
    If mode == 'center', centers the slab with vacuum split top and bottom.
    """
    lat = np.array(structure.lattice.matrix, dtype=float)
    if np.max(np.abs(lat[:2, 2])) > 1e-6 or np.max(np.abs(lat[2, :2])) > 1e-6:
        raise RuntimeError("Slab must be oriented with a,b in xy and c along z")
    z = np.asarray(structure.cart_coords[:, 2], dtype=float)
    thickness = float(np.ptp(z))
    new_c = thickness + vacuum_angstrom
    lat[2] = [0.0, 0.0, new_c]
    f = np.array(structure.frac_coords, dtype=float)
    f[:, :2] %= 1.0
    if mode == "top":
        f[:, 2] = np.clip((z - z.min()) / new_c, 0.0, 1.0)
    elif mode == "center":
        f[:, 2] = (z - z.min() + vacuum_angstrom / 2.0) / new_c
    else:
        raise ValueError(f"Unknown vacuum mode: {mode}")
    return Structure(Lattice(lat), structure.species, f)


def calculate_checks(
    structure: Structure,
    plane_tolerance: float,
    vacuum_target: float,
    freeze_bottom_planes: int,
    overlap_threshold: float,
    vacuum_mode: str = "top",
) -> tuple[list[list[bool]], dict, list[str]]:
    """Create selective-dynamics flags and validate the produced slab."""
    z = np.asarray(structure.cart_coords[:, 2], dtype=float)
    planes = group_atomic_planes(z, plane_tolerance)
    if len(planes) <= freeze_bottom_planes:
        raise RuntimeError(
            f"Only {len(planes)} distinct atomic z groups; cannot fix "
            f"{freeze_bottom_planes} and leave surface atoms mobile"
        )

    frozen = set(i for plane in planes[:freeze_bottom_planes] for i in plane)
    flags = [[i not in frozen] * 3 for i in range(len(structure))]
    free = sorted(set(range(len(structure))) - frozen)
    frozen_sorted = sorted(frozen)

    lattice = structure.lattice
    matrix = np.asarray(lattice.matrix)
    thickness = float(np.ptp(z))
    vacuum = float(lattice.c - thickness)
    vacuum_bottom = float(np.min(z))
    vacuum_top = float(lattice.c - np.max(z))
    if abs(vacuum - vacuum_target) > 1e-5:
        raise RuntimeError(f"Vacuum gap {vacuum:.6f} A disagrees with {vacuum_target:.6f} A")
    if np.linalg.det(matrix) <= 0:
        raise RuntimeError("Output cell is left-handed")
    if not np.all((structure.frac_coords >= -1e-6) & (structure.frac_coords <= 1 + 1e-6)):
        raise RuntimeError("Atomic coordinates fall outside the unit cell")

    n = len(structure)
    if n > 1:
        distances = np.asarray(structure.distance_matrix)
        iu = np.triu_indices(n, 1)
        min_dist = float(np.min(distances[iu]))
    else:
        min_dist = float('inf')
    if min_dist < overlap_threshold:
        raise RuntimeError(
            f"Detected unphysically close periodic atoms: {min_dist:.4f} A "
            f"(< {overlap_threshold} A)"
        )

    counts = species_counts(structure)
    # Nominal Li(+1), Co(+3), O(-2) bookkeeping is a termination screen,
    # not a calculated Bader charge or prediction of electron distribution.
    nominal_ionic_charge = (
        counts.get("Li", 0) + 3 * counts.get("Co", 0)
        - 2 * counts.get("O", 0)
    )
    warnings = []
    if nominal_ionic_charge != 0:
        warnings.append(
            f"Nominal Li(+1)/Co(+3)/O(-2) ionic charge sum is "
            f"{nominal_ionic_charge:+d} elementary charges per slab cell. "
            "Check polar/nonstoichiometric termination and compensation before DFT."
        )
    if not (counts.get("Li", 0) == counts.get("Co", 0)
            and counts.get("O", 0) == 2 * counts.get("Co", 0)):
        warnings.append(
            "Slab is nonstoichiometric relative to LiCoO2. This can be valid "
            "for a particular termination, but check surface chemistry, charge, "
            "dipole and reference-energy methodology."
        )
    warnings.extend([
        "'3 planes' is a nominal pymatgen Miller-plane thickness request, "
        "NOT proof of three complete LiCoO2 chemical layers.",
        "Atomic z groups use a tolerance-based geometric clustering; "
        "they are not automatically chemically complete layers.",
        "Inspect both surface terminations visually before DFT. For an "
        "asymmetric slab, evaluate dipole correction (IDIPOL=3/LDIPOL) "
        "and vacuum convergence.",
        f"Fixing bottom {freeze_bottom_planes} z groups is a starting constraint choice. "
        "Verify that fixed atoms are bulk-like and enough atoms remain free.",
    ])
    if min_dist < 1.2:
        warnings.append("Very short interatomic separation: inspect geometry for defects.")
    if n > 400:
        warnings.append("Large slab: benchmark 40-core VASP memory and MPI efficiency first.")

    plane_info = []
    for k, indices in enumerate(planes, start=1):
        plane_info.append({
            "group_from_bottom": k,
            "mean_z_angstrom": float(np.mean(z[indices])),
            "atom_count": len(indices),
            "composition": species_counts(structure, indices),
            "poscar_atom_indices_1_based": [i + 1 for i in sorted(indices)],
            "frozen": k <= freeze_bottom_planes,
        })
    report = {
        "atom_count": n,
        "species_order": list(SPECIES_ORDER),
        "composition": counts,
        "nominal_ionic_charge_e": nominal_ionic_charge,
        "lattice_abc_angstrom": list(map(float, lattice.abc)),
        "lattice_angles_degrees": list(map(float, lattice.angles)),
        "lattice_vectors_angstrom": matrix.tolist(),
        "slab_atomic_span_angstrom": thickness,
        "vacuum_geometric_gap_angstrom": vacuum,
        "vacuum_mode": vacuum_mode,
        "vacuum_at_top_angstrom": vacuum_top,
        "vacuum_at_bottom_angstrom": vacuum_bottom,
        "minimum_periodic_atom_distance_angstrom": min_dist,
        "geometric_atomic_z_group_count": len(planes),
        "atomic_z_groups": plane_info,
        "frozen_atom_count": len(frozen_sorted),
        "free_atom_count": len(free),
        "frozen_atom_indices_1_based": [i + 1 for i in frozen_sorted],
        "frozen_atom_species": species_counts(structure, frozen_sorted),
        "free_atom_species": species_counts(structure, free),
        "bottommost_group_species": plane_info[0]["composition"],
        "topmost_group_species": plane_info[-1]["composition"],
        "dipol_center_fractional_z_approx": float(np.mean(z) / lattice.c),
    }
    return flags, report, warnings


def write_qe_blocks(path: Path, structure: Structure, flags: list[list[bool]]) -> None:
    """Optional legacy QE coordinate blocks; NOT a complete pw.x input."""
    lines = ["CELL_PARAMETERS (angstrom)"]
    for vec in structure.lattice.matrix:
        lines.append("  " + "  ".join(f"{x:17.10f}" for x in vec))
    lines.append("\nATOMIC_POSITIONS (angstrom)")
    for site, flag in zip(structure, flags):
        coords = " ".join(f"{x:16.10f}" for x in site.coords)
        qe_flags = " ".join("1" if allowed else "0" for allowed in flag)
        lines.append(f"{site.specie.symbol:<2}  {coords}   {qe_flags}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build a VASP LiCoO2 (104) slab from optimized bulk CONTCAR."
    )
    p.add_argument("--bulk", type=Path, default=Path(DEFAULT_BULK),
                   help=f"Relaxed bulk CONTCAR or CIF (default: {DEFAULT_BULK})")
    p.add_argument("--output-dir", type=Path, default=Path("./lco_104_slab"))
    p.add_argument("--miller", type=int, nargs=3, default=(1, 0, 4), metavar=("H", "K", "L"))
    p.add_argument("--thickness-mode", choices=("planes", "angstrom"), default="planes",
                   help="'planes' means nominal Miller-plane spacings, not chemical layers")
    p.add_argument("--slab-size", type=float, default=3.0,
                   help="Minimum thickness in requested mode (default: 3 nominal planes)")
    p.add_argument("--vacuum", type=float, default=18.0,
                   help="Geometric gap between periodic slab images in Angstrom (default: 18)")
    p.add_argument("--vacuum-mode", choices=("top", "center"), default="top",
                   help="Where to place the vacuum gap: 'top' puts all vacuum above the slab (bottom layer at z=0); 'center' splits vacuum equally top and bottom (default: top)")
    p.add_argument("--supercell", nargs=3, type=int, default=(3, 2, 1),
                   metavar=("A", "B", "C"), help="In-plane expansion; C must be 1")
    p.add_argument("--termination-index", type=int, default=0,
                   help="Choose from printed termination inventory (default: 0)")
    p.add_argument("--list-terminations", action="store_true",
                   help="Only print termination inventory; do not generate files")
    p.add_argument("--freeze-bottom-planes", type=int, default=2,
                   help="Freeze bottom N *atomic z clusters* (default: 2)")
    p.add_argument("--plane-tolerance", type=float, default=0.35,
                   help="Tolerance for grouping atomic z coordinates in Angstrom")
    p.add_argument("--minimum-distance", type=float, default=0.8,
                   help="Fail if atoms are closer than this periodic distance (Angstrom)")
    p.add_argument("--keep-original-origin", action="store_true",
                   help="Do not translate the lowest-z Li onto fractional x=y=0")
    p.add_argument("--skip-reduction", action="store_true",
                   help="Skip optional unimodular Gauss reduction of surface basis")
    p.add_argument("--qe-blocks", action="store_true",
                   help="Additionally write QE_blocks.in with 0/1 frozen/free flags")
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if not args.bulk.is_file():
        raise FileNotFoundError(f"Bulk structure not found: {args.bulk}")
    if args.slab_size <= 0 or args.vacuum <= 0 or args.plane_tolerance <= 0:
        raise ValueError("Slab size, vacuum and plane tolerance must be positive")
    if args.minimum_distance <= 0 or args.freeze_bottom_planes < 0:
        raise ValueError("Minimum distance must be positive; freeze count >= 0")
    if any(x <= 0 for x in args.supercell) or args.supercell[2] != 1:
        raise ValueError("supercell must be positive and third multiplier must equal 1")
    if args.miller == [0, 0, 0]:
        raise ValueError("Miller indices cannot all be zero")

    bulk = Structure.from_file(str(args.bulk))
    if set(site.specie.symbol for site in bulk) != set(SPECIES_ORDER):
        raise ValueError("Expected bulk containing exactly Li, Co and O")
    if np.linalg.det(bulk.lattice.matrix) <= 0:
        raise ValueError("Bulk lattice must be right-handed")

    hkl = tuple(args.miller)
    dhkl = float(bulk.lattice.d_hkl(hkl))
    # In unit-plane mode, pymatgen interprets *both* slab and vacuum
    # minima in plane units. Convert the user's requested vacuum (Angstrom)
    # to plane units, then adjust the final geometric gap precisely below.
    plane_mode = args.thickness_mode == "planes"
    vacuum_param = float(math.ceil(args.vacuum / dhkl)) if plane_mode else args.vacuum
    generator = SlabGenerator(
        initial_structure=bulk,
        miller_index=hkl,
        min_slab_size=args.slab_size,
        min_vacuum_size=vacuum_param,
        center_slab=True,
        in_unit_planes=plane_mode,
        primitive=True,
        reorient_lattice=True,
    )
    slabs = generator.get_slabs()
    if not slabs:
        raise RuntimeError(f"No slab terminations returned for {hkl}")

    inventory = termination_inventory(slabs, args.plane_tolerance)
    print(f"Bulk: {args.bulk.resolve()} ({len(bulk)} atoms)")
    print(f"Miller index: {hkl}; d_hkl={dhkl:.4f} A")
    print(f"Thickness: minimum {args.slab_size:g} {args.thickness_mode}")
    print("Available terminations (before supercell):")
    for item in inventory:
        print(f"  [{item['termination_index']}] shift={item['pymatgen_shift']:.5f} "
              f"atoms={item['atoms_in_unexpanded_cell']} "
              f"bottom={item['bottom_plane_species']} "
              f"top={item['top_plane_species']} "
              f"total={item['composition']}")
    if args.list_terminations:
        return 0
    if not 0 <= args.termination_index < len(slabs):
        raise ValueError(f"Termination index must be in [0, {len(slabs) - 1}]")

    chosen = slabs[args.termination_index].get_orthogonal_c_slab()
    if not args.skip_reduction:
        chosen = compact_in_plane(chosen)
    chosen = orient_surface_xy(chosen)
    chosen.make_supercell(list(args.supercell))
    if not args.keep_original_origin:
        chosen, li_origin_index = translate_li_xy_origin(chosen)
    else:
        li_origin_index = None
    final = recenter_with_vacuum(chosen, args.vacuum, mode=args.vacuum_mode)
    order = {symbol: i for i, symbol in enumerate(SPECIES_ORDER)}
    final = final.get_sorted_structure(key=lambda site: order[site.specie.symbol])

    flags, report, warnings = calculate_checks(
        final, args.plane_tolerance, args.vacuum,
        args.freeze_bottom_planes, args.minimum_distance,
        vacuum_mode=args.vacuum_mode,
    )
    try:
        pymatgen_version = version("pymatgen")
    except PackageNotFoundError:
        pymatgen_version = "unknown"
    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pymatgen_version": pymatgen_version,
        "source_bulk_file": str(args.bulk.resolve()),
        "source_bulk_lattice_abc_angstrom": list(map(float, bulk.lattice.abc)),
        "miller_index": list(hkl),
        "miller_spacing_angstrom": dhkl,
        "thickness_request": {"mode": args.thickness_mode, "value": args.slab_size},
        "vacuum_requested_angstrom": args.vacuum,
        "vacuum_mode": args.vacuum_mode,
        "supercell": list(args.supercell),
        "termination_index": args.termination_index,
        "selected_termination": inventory[args.termination_index],
        "all_terminations": inventory,
        "li_origin_pre_sort_index_1_based": (
            li_origin_index + 1 if li_origin_index is not None else None
        ),
        "fixed_bottom_z_groups_requested": args.freeze_bottom_planes,
        "plane_tolerance_angstrom": args.plane_tolerance,
        "geometry": report,
        "warnings": warnings,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    poscar_path = args.output_dir / "POSCAR"
    Poscar(final, comment="LiCoO2 (104) candidate; verify termination",
           selective_dynamics=flags, sort_structure=False).write_file(poscar_path)
    CifWriter(final, refine_struct=False).write_file(str(args.output_dir / "slab_104.cif"))
    (args.output_dir / "slab_report.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    text_lines = [
        "LiCoO2 (104) SLAB GEOMETRY VALIDATION (PRE-RELAXATION)",
        f"Source bulk: {metadata['source_bulk_file']}",
        f"Miller index: {hkl}",
        f"Requested: {args.slab_size:g} {args.thickness_mode}; "
        f"{args.vacuum:g} A vacuum; supercell {tuple(args.supercell)}",
        f"Selected termination: {args.termination_index} (not necessarily lowest-energy)",
        f"Species/POTCAR order: {' '.join(SPECIES_ORDER)}",
        f"Composition: {report['composition']}; atoms: {report['atom_count']}",
        f"Nominal ionic charge proxy [e]: {report['nominal_ionic_charge_e']:+d}",
        f"a,b,c [A]: {report['lattice_abc_angstrom']}",
        f"alpha,beta,gamma [degrees]: {report['lattice_angles_degrees']}",
        f"Atomic thickness [A]: {report['slab_atomic_span_angstrom']:.4f}",
        f"Geometric vacuum gap [A]: {report['vacuum_geometric_gap_angstrom']:.4f} (mode: {args.vacuum_mode})",
        f"  - Vacuum above top surface [A]: {report['vacuum_at_top_angstrom']:.4f}",
        f"  - Vacuum below bottom surface [A]: {report['vacuum_at_bottom_angstrom']:.4f}",
        f"Minimum pair distance [A]: {report['minimum_periodic_atom_distance_angstrom']:.4f}",
        f"Atomic z groups: {report['geometric_atomic_z_group_count']}",
        f"Frozen: {report['frozen_atom_count']}; free: {report['free_atom_count']}",
        f"Frozen atom indices in POSCAR (1-based): {report['frozen_atom_indices_1_based']}",
        f"Frozen species: {report['frozen_atom_species']}",
        f"Free species: {report['free_atom_species']}",
        "\nAtomic z groups from bottom (not necessarily full chemical layers):",
    ]
    for item in report["atomic_z_groups"]:
        text_lines.append(
            f"  {item['group_from_bottom']:>2}: z={item['mean_z_angstrom']:8.4f} A "
            f"atoms={item['atom_count']:3} comp={item['composition']} "
            f"fixed={item['frozen']} indices={item['poscar_atom_indices_1_based']}"
        )
    text_lines.append("\nWarnings and required review:")
    for warning in warnings:
        text_lines.append(f"  - {warning}")
    text_lines.extend([
        "\nNext steps:",
        "  1. Inspect POSCAR/CIF in VESTA; verify (104), chemistry and layer count.",
        "  2. Confirm fixed indices reflect physically sensible bulk-like layers.",
        "  3. Check stoichiometry, dipole and termination suitability.",
        "  4. Supply licensed PBE64 Li_sv/Co_pv/O POTCAR in that species order.",
        "  5. Use separate relaxation INCAR (ISIF=2), converged k-grid, and SLURM.",
        "  6. Converge thickness, vacuum, k-mesh and ENCUT before final adsorption energies.",
    ])
    (args.output_dir / "slab_report.txt").write_text(
        "\n".join(text_lines) + "\n", encoding="utf-8"
    )
    if args.qe_blocks:
        write_qe_blocks(args.output_dir / "QE_blocks.in", final, flags)

    print("\n".join(text_lines[:17]))
    print("\nGenerated files:")
    for name in ["POSCAR", "slab_104.cif", "slab_report.txt", "slab_report.json"]:
        print(" ", (args.output_dir / name).resolve())
    if args.qe_blocks:
        print(" ", (args.output_dir / "QE_blocks.in").resolve())
    print("Review slab_report.txt and inspect structure before VASP submission.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
