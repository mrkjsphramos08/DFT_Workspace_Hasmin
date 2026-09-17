"""Build adsorbed M4 cluster structures on an LCO slab.

Places each relaxed tetrahedral cluster (oriented base-down) at the chosen
adsorption sites, 2.0 A above the topmost surface atoms, and writes
``adsorbed_<cluster-stem>_<site>.cif`` under the slab profile's namespace.
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

METALS_DIR = PROJECT_ROOT / "05_metals"
ADSORPTION_HEIGHT = 2.0  # Å distance between the top surface layer and cluster base
SURFACE_Z_TOLERANCE = 0.8
PAIR_MAX_DISTANCE = 4.5


def sanitize_atoms(atoms_obj):
    """Drop leftover CIF metadata so ASE can write the structure cleanly."""
    return Atoms(
        symbols=atoms_obj.get_chemical_symbols(),
        positions=atoms_obj.get_positions(),
        cell=atoms_obj.get_cell(),
        pbc=atoms_obj.get_pbc(),
    )


def orient_tetrahedron_base_down(cluster):
    """Rotate a 4-atom tetrahedron so its triangular base points downward."""
    pos = cluster.get_positions()
    if len(pos) != 4:
        return cluster

    best_angle = float("inf")
    best_normal = None
    indices = [0, 1, 2, 3]

    for apex_idx in indices:
        base_indices = [j for j in indices if j != apex_idx]
        a, b, c = pos[base_indices]
        apex = pos[apex_idx]

        normal = np.cross(b - a, c - a)
        if np.dot(normal, apex - a) < 0:
            normal = -normal

        norm_len = np.linalg.norm(normal)
        if norm_len > 1e-5:
            normal /= norm_len

        angle = np.arccos(np.clip(np.dot(normal, [0.0, 0.0, 1.0]), -1.0, 1.0))
        if angle < best_angle:
            best_angle = angle
            best_normal = normal

    if best_normal is not None:
        cluster.rotate(best_normal, (0, 0, 1), rotate_cell=False)

    return cluster


def get_top_surface_indices(slab, element, z_tol=SURFACE_Z_TOLERANCE):
    """Return the top-most atoms of a chosen element in the surface layer."""
    valid_mask = slab.get_scaled_positions()[:, 2] < 0.8
    valid_indices = np.where(valid_mask)[0]
    matching = [i for i in valid_indices if slab.symbols[i] == element]
    if not matching:
        return []

    max_z = np.max(slab.positions[matching, 2])
    return [i for i in matching if slab.positions[i, 2] > (max_z - z_tol)]


def find_best_pair_midpoint(slab, site_indices, center_xy, max_pair_dist=PAIR_MAX_DISTANCE):
    """Pick the nearest-neighbor pair midpoint closest to the slab center."""
    best_midpoint = None
    min_dist = float("inf")

    for i in range(len(site_indices)):
        for j in range(i + 1, len(site_indices)):
            idx1, idx2 = site_indices[i], site_indices[j]
            pair_distance = slab.get_distance(idx1, idx2, mic=True)
            if pair_distance >= max_pair_dist:
                continue

            vec_12 = slab.get_distance(idx1, idx2, vector=True, mic=True)
            midpoint_xy = (slab.positions[idx1] + vec_12 / 2.0)[:2]
            dist_to_center = np.linalg.norm(midpoint_xy - center_xy)

            if dist_to_center < min_dist:
                min_dist = dist_to_center
                best_midpoint = midpoint_xy

    return best_midpoint


def build_target_sites(slab):
    """Assemble the final list of adsorption targets used for this slab."""
    scaled_pos = slab.get_scaled_positions()
    valid_mask = scaled_pos[:, 2] < 0.8
    surface_z = np.max(slab.positions[valid_mask, 2])

    cell = slab.get_cell()
    slab_xy_center = (cell[0][:2] + cell[1][:2]) / 2.0

    li_indices = get_top_surface_indices(slab, "Li")
    co_indices = get_top_surface_indices(slab, "Co")
    o_indices = get_top_surface_indices(slab, "O")

    missing = []
    if not li_indices:
        missing.append("Li")
    if not co_indices:
        missing.append("Co")
    if not o_indices:
        missing.append("O")
    if missing:
        raise ValueError(f"Could not find surface {', '.join(missing)} atoms in the slab.")

    best_li_idx = min(li_indices, key=lambda i: np.linalg.norm(slab.positions[i][:2] - slab_xy_center))
    best_co_idx = min(co_indices, key=lambda i: np.linalg.norm(slab.positions[i][:2] - slab_xy_center))
    best_o_idx = min(o_indices, key=lambda i: np.linalg.norm(slab.positions[i][:2] - slab_xy_center))

    best_lili_midpoint = find_best_pair_midpoint(slab, li_indices, slab_xy_center)
    best_coco_midpoint = find_best_pair_midpoint(slab, co_indices, slab_xy_center)

    best_lico_midpoint = None
    min_midpoint_dist = float("inf")
    for li_i in li_indices:
        for co_j in co_indices:
            vec_ij = slab.get_distance(li_i, co_j, vector=True, mic=True)
            midpoint_xy = (slab.positions[li_i] + vec_ij / 2.0)[:2]
            dist_to_center = np.linalg.norm(midpoint_xy - slab_xy_center)
            if dist_to_center < min_midpoint_dist:
                min_midpoint_dist = dist_to_center
                best_lico_midpoint = midpoint_xy

    target_sites = {
        "li": slab.positions[best_li_idx][:2],
        "co": slab.positions[best_co_idx][:2],
        "lico": best_lico_midpoint,
        "centered": slab_xy_center,
        "oxygen": slab.positions[best_o_idx][:2],
        "lili": best_lili_midpoint,
        "coco": best_coco_midpoint,
    }

    ordered_target_sites = {site: target_sites.get(site) for site in SITE_ORDER if site in target_sites}
    return surface_z, ordered_target_sites


def place_cluster_on_site(cluster_path, slab, surface_z, site_name, target_xy, slab_name):
    """Build one adsorbed structure by placing the cluster at the chosen site."""
    cluster = read(cluster_path)
    cluster = orient_tetrahedron_base_down(cluster)

    cluster_z_min = np.min(cluster.positions[:, 2])
    base_mask = cluster.positions[:, 2] < (cluster_z_min + 0.5)
    base_xy_center = np.mean(cluster.positions[base_mask, :2], axis=0)

    xy_shift = target_xy - base_xy_center
    z_shift = (surface_z + ADSORPTION_HEIGHT) - cluster_z_min
    cluster.positions += np.array([xy_shift[0], xy_shift[1], z_shift])

    adsorbed = slab.copy()
    adsorbed.extend(cluster)
    clean_adsorbed = sanitize_atoms(adsorbed)

    site_dir = slab_output_root(cluster_path.parent, slab_name) / site_name
    site_dir.mkdir(parents=True, exist_ok=True)

    out_path = site_dir / f"adsorbed_{cluster_path.stem}_{site_name}.cif"
    write(out_path, clean_adsorbed)
    print(f"Saved [{site_name.upper()}]: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate tetrahedral clusters adsorbed on an LCO slab.")
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

    tetra_files = [
        path for path in METALS_DIR.rglob("*_latest.cif")
        if "Tetra" in path.parts and "Adsorbed" not in path.parts
    ]
    if not tetra_files:
        print("No tetrahedral cluster CIF files (*_latest.cif) found in 05_metals/")
        return

    print(f"Found {len(tetra_files)} tetrahedral cluster(s) to process.\n")

    for cluster_path in tetra_files:
        for site_name, target_xy in target_sites.items():
            if target_xy is None:
                continue
            place_cluster_on_site(cluster_path, slab, surface_z, site_name, target_xy, args.slab)


if __name__ == "__main__":
    main()