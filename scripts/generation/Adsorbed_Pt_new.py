"""Build adsorbed single-atom Pt structures on an LCO slab.

Places the Pt atom from the Pt1 folder at each site defined in
``adsorption_sites.SITE_SPECS``, 2.5 A above the topmost surface atoms, and writes
``adsorbed_<pt1-stem>_<site>.cif`` under the slab profile's namespace.
Surface-site selection works in fractional ``z`` because slab CIFs wrap across
the periodic boundary.
"""
import argparse
from itertools import combinations, product
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read, write

from adsorption_config import DEFAULT_SLAB, SLAB_CHOICES, get_slab_profile, slab_output_root
from adsorption_sites import SITE_SPECS

SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = next(
    (p for p in [SCRIPT_PATH.parent] + list(SCRIPT_PATH.parents) if (p / "05_metals").exists()),
    SCRIPT_PATH.parent,
)

PT1_DIR = PROJECT_ROOT / "05_metals" / "Pt1"
ADSORPTION_HEIGHT = 2.5  # Å between the top surface layer and the Pt atom
SURFACE_Z_TOLERANCE = 0.8
SURFACE_FRAC_Z_MAX = 0.8  # ignore atoms wrapped across the periodic boundary


def sanitize_atoms(atoms_obj):
    """Drop leftover CIF metadata so ASE can write the structure cleanly."""
    return Atoms(
        symbols=atoms_obj.get_chemical_symbols(),
        positions=atoms_obj.get_positions(),
        cell=atoms_obj.get_cell(),
        pbc=atoms_obj.get_pbc(),
    )


def get_top_surface_indices(slab, element, z_tol=SURFACE_Z_TOLERANCE):
    """Return the top-most atoms of a chosen element in the surface layer."""
    valid_mask = slab.get_scaled_positions()[:, 2] < SURFACE_FRAC_Z_MAX
    matching = [i for i in np.where(valid_mask)[0] if slab.symbols[i] == element]
    if not matching:
        return []
    max_z = np.max(slab.positions[matching, 2])
    return [i for i in matching if slab.positions[i, 2] > (max_z - z_tol)]


def closest_to_center(slab, indices, center_xy):
    """Return the atom index whose xy position is nearest the slab center."""
    return min(indices, key=lambda i: np.linalg.norm(slab.positions[i][:2] - center_xy))


def surface_pairs(slab, surface, elements, max_dist):
    """Pairs of surface atoms for a bridge site, optionally limited by distance."""
    el_a, el_b = elements
    if el_a == el_b:
        pairs = combinations(surface[el_a], 2)
    else:
        pairs = product(surface[el_a], surface[el_b])
    return [
        (i, j) for i, j in pairs
        if max_dist is None or slab.get_distance(i, j, mic=True) < max_dist
    ]


def find_best_midpoint(slab, pairs, center_xy):
    """Pick the (i, j) pair midpoint closest to the slab center (xy), or None."""
    best_midpoint, min_dist = None, float("inf")
    for i, j in pairs:
        vec = slab.get_distance(i, j, vector=True, mic=True)
        midpoint_xy = (slab.positions[i] + vec / 2.0)[:2]
        dist = np.linalg.norm(midpoint_xy - center_xy)
        if dist < min_dist:
            min_dist, best_midpoint = dist, midpoint_xy
    return best_midpoint


def resolve_site(slab, spec, surface, center_xy):
    """Return the xy target for one site spec (None if it cannot be found)."""
    kind = spec["kind"]
    if kind == "center":
        return center_xy
    if kind == "atop":
        idx = closest_to_center(slab, surface[spec["elements"][0]], center_xy)
        return slab.positions[idx][:2]
    if kind == "bridge":
        pairs = surface_pairs(slab, surface, spec["elements"], spec.get("max_dist"))
        return find_best_midpoint(slab, pairs, center_xy)
    raise ValueError(f"Unknown site kind: {kind!r}")


def build_target_sites(slab):
    """Assemble the adsorption targets for this slab, in SITE_SPECS order."""
    valid_mask = slab.get_scaled_positions()[:, 2] < SURFACE_FRAC_Z_MAX
    surface_z = np.max(slab.positions[valid_mask, 2])

    cell = slab.get_cell()
    center_xy = (cell[0][:2] + cell[1][:2]) / 2.0

    needed = sorted({el for spec in SITE_SPECS.values() for el in spec["elements"]})
    surface = {el: get_top_surface_indices(slab, el) for el in needed}
    missing = [el for el, idx in surface.items() if not idx]
    if missing:
        raise ValueError(f"Could not find surface {', '.join(missing)} atoms in the slab.")

    sites = {}
    for name, spec in SITE_SPECS.items():
        xy = resolve_site(slab, spec, surface, center_xy)
        if xy is None:
            print(f"WARNING: no valid pair found for site '{name}' "
                  f"(max_dist={spec.get('max_dist')}); skipping.")
            continue
        sites[name] = xy
    return surface_z, sites


def place_atom_on_site(pt_path, slab, surface_z, site_name, target_xy, slab_name):
    """Build one adsorbed structure by placing the Pt atom at the chosen site."""
    pt = read(pt_path)
    if len(pt) != 1:
        raise ValueError(f"{pt_path} contains {len(pt)} atoms; expected a single Pt atom.")

    pt.positions[0] = [target_xy[0], target_xy[1], surface_z + ADSORPTION_HEIGHT]

    adsorbed = slab.copy()
    adsorbed.extend(pt)

    site_dir = slab_output_root(pt_path.parent, slab_name) / site_name
    site_dir.mkdir(parents=True, exist_ok=True)

    out_path = site_dir / f"adsorbed_{pt_path.stem}_{site_name}.cif"
    write(out_path, sanitize_atoms(adsorbed))
    print(f"Saved [{site_name.upper()}]: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate a single Pt atom adsorbed on an LCO slab.")
    parser.add_argument("--slab", choices=SLAB_CHOICES, default=DEFAULT_SLAB,
                        help=f"Slab profile to use (default: {DEFAULT_SLAB})")
    parser.add_argument("--sites", nargs="+", choices=list(SITE_SPECS),
                        help="Only generate these sites (default: all)")
    args = parser.parse_args()

    slab_path = get_slab_profile(args.slab)["slab_path"]
    if not slab_path.exists():
        raise FileNotFoundError(f"Slab file missing for profile '{args.slab}': {slab_path}")

    slab = read(slab_path)
    surface_z, target_sites = build_target_sites(slab)
    if args.sites:
        target_sites = {s: xy for s, xy in target_sites.items() if s in args.sites}

    pt_files = [p for p in PT1_DIR.rglob("*_latest.cif") if "Adsorbed" not in p.parts]
    if not pt_files:
        print(f"No Pt1 CIF files (*_latest.cif) found in {PT1_DIR}")
        return

    print(f"Found {len(pt_files)} Pt1 structure(s) to process.\n")
    for pt_path in pt_files:
        for site_name, target_xy in target_sites.items():
            place_atom_on_site(pt_path, slab, surface_z, site_name, target_xy, args.slab)


if __name__ == "__main__":
    main()