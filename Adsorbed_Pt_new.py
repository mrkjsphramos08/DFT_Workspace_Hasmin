"""Build adsorbed single-atom Pt structures on an LCO slab.

Places the Pt atom from the Pt1 folder at each chosen adsorption site,
2.0 A above the topmost surface atoms, and writes
``adsorbed_<pt1-stem>_<site>.cif`` under the slab profile's namespace.
Surface-site selection works in fractional ``z`` because slab CIFs wrap across
the periodic boundary.
"""
import argparse
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read, write

from adsorption_config import DEFAULT_SLAB, SLAB_CHOICES, get_slab_profile, slab_output_root
from adsorption_sites import SITE_ORDER

# -------------------------------------------------------------------------
# Configuration and path resolution
# -------------------------------------------------------------------------
SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = next(
    (p for p in [SCRIPT_PATH.parent] + list(SCRIPT_PATH.parents) if (p / "05_metals").exists()),
    SCRIPT_PATH.parent,
)

PT1_DIR = PROJECT_ROOT / "05_metals" / "Pt1"
ADSORPTION_HEIGHT = 2.0  # Å between the top surface layer and the Pt atom
SURFACE_Z_TOLERANCE = 0.8
SURFACE_FRAC_Z_MAX = 0.8  # ignore atoms wrapped across the periodic boundary
PAIR_MAX_DISTANCE = 4.5


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


def find_best_midpoint(slab, pairs, center_xy):
    """Pick the (i, j) pair midpoint closest to the slab center."""
    best_midpoint, min_dist = None, float("inf")
    for i, j in pairs:
        vec = slab.get_distance(i, j, vector=True, mic=True)
        midpoint_xy = (slab.positions[i] + vec / 2.0)[:2]
        dist = np.linalg.norm(midpoint_xy - center_xy)
        if dist < min_dist:
            min_dist, best_midpoint = dist, midpoint_xy
    return best_midpoint


def same_element_pairs(slab, indices, max_dist=PAIR_MAX_DISTANCE):
    """Nearest-neighbor pairs among the given indices."""
    return [
        (i, j)
        for n, i in enumerate(indices)
        for j in indices[n + 1:]
        if slab.get_distance(i, j, mic=True) < max_dist
    ]


def build_target_sites(slab):
    """Assemble the final list of adsorption targets used for this slab."""
    valid_mask = slab.get_scaled_positions()[:, 2] < SURFACE_FRAC_Z_MAX
    surface_z = np.max(slab.positions[valid_mask, 2])

    cell = slab.get_cell()
    center_xy = (cell[0][:2] + cell[1][:2]) / 2.0

    surface = {el: get_top_surface_indices(slab, el) for el in ("Li", "Co", "O")}
    missing = [el for el, idx in surface.items() if not idx]
    if missing:
        raise ValueError(f"Could not find surface {', '.join(missing)} atoms in the slab.")

    lico_pairs = [(li, co) for li in surface["Li"] for co in surface["Co"]]

    target_sites = {
        "li": slab.positions[closest_to_center(slab, surface["Li"], center_xy)][:2],
        "co": slab.positions[closest_to_center(slab, surface["Co"], center_xy)][:2],
        "lico": find_best_midpoint(slab, lico_pairs, center_xy),
        "centered": center_xy,
        "oxygen": slab.positions[closest_to_center(slab, surface["O"], center_xy)][:2],
        "lili": find_best_midpoint(slab, same_element_pairs(slab, surface["Li"]), center_xy),
        "coco": find_best_midpoint(slab, same_element_pairs(slab, surface["Co"]), center_xy),
    }

    ordered = {site: target_sites[site] for site in SITE_ORDER if site in target_sites}
    return surface_z, ordered


def place_atom_on_site(pt_path, slab, surface_z, site_name, target_xy, slab_name):
    """Build one adsorbed structure by placing the Pt atom at the chosen site."""
    pt = read(pt_path)
    if len(pt) != 1:
        raise ValueError(f"{pt_path} contains {len(pt)} atoms; expected a single Pt atom.")

    pt.positions[0] = [target_xy[0], target_xy[1], surface_z + ADSORPTION_HEIGHT]

    adsorbed = slab.copy()
    adsorbed.extend(pt)
    clean_adsorbed = sanitize_atoms(adsorbed)

    site_dir = slab_output_root(pt_path.parent, slab_name) / site_name
    site_dir.mkdir(parents=True, exist_ok=True)

    out_path = site_dir / f"adsorbed_{pt_path.stem}_{site_name}.cif"
    write(out_path, clean_adsorbed)
    print(f"Saved [{site_name.upper()}]: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate a single Pt atom adsorbed on an LCO slab.")
    parser.add_argument(
        "--slab",
        choices=SLAB_CHOICES,
        default=DEFAULT_SLAB,
        help=f"Slab profile to use (default: {DEFAULT_SLAB})",
    )
    args = parser.parse_args()

    slab_path = get_slab_profile(args.slab)["slab_path"]
    if not slab_path.exists():
        raise FileNotFoundError(f"Slab file missing for profile '{args.slab}': {slab_path}")

    slab = read(slab_path)
    surface_z, target_sites = build_target_sites(slab)

    pt_files = [p for p in PT1_DIR.rglob("*_latest.cif") if "Adsorbed" not in p.parts]
    if not pt_files:
        print(f"No Pt1 CIF files (*_latest.cif) found in {PT1_DIR}")
        return

    print(f"Found {len(pt_files)} Pt1 structure(s) to process.\n")

    for pt_path in pt_files:
        for site_name, target_xy in target_sites.items():
            if target_xy is None:
                continue
            place_atom_on_site(pt_path, slab, surface_z, site_name, target_xy, args.slab)


if __name__ == "__main__":
    main()