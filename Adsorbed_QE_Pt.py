"""Generate pw.x inputs for the adsorbed single-atom Pt systems.

Uses the bare-atom input in the ``Pt1`` directory as a template,
rewrites ``nat``/``ntyp`` and the ATOMIC_SPECIES block, freezes the bottom half
of the slab with ``0 0 0`` flags, and appends the Hubbard and k-point cards.
Sites that already contain a QE output are skipped, so re-running is idempotent.
"""
import argparse
import re
from pathlib import Path
from ase.io import read

from adsorption_config import DEFAULT_SLAB, SLAB_CHOICES, get_slab_profile
from adsorption_sites import SITE_ORDER, VALID_SITES
from lco.paths import PROJECT_ROOT
from lco.qe_input import (
    PSEUDO_DB,
    CARDS_TO_STRIP,
    HUBBARD_BLOCK,
    build_atomic_species_block,
    strip_old_cards,
    make_prefix_directive,
    SURFACE_KPOINTS,
    variant_prefix,
)

# -------------------------------------------------------------------------
# Configuration & Path Resolution
# -------------------------------------------------------------------------
SCRIPT_PATH = Path(__file__).resolve()
PT1_DIR = PROJECT_ROOT / "05_metals" / "Pt1"

NUM_ADSORBATE_ATOMS = 1  # single Pt atom, always the last atom in the CIF
Z_CUTOFF_FRACTION = 0.5  # Freeze bottom 50% of slab height

# Additional metadata needed for unique prefixes
SYSTEM_TAG = "Pt1"


def find_templates_in(pt1_dir):
    """Find all bare-atom template .in files directly in the Pt1 directory."""
    if not pt1_dir.exists():
        return []

    # Ignore any 'adsorbed' input files; sort so the choice is deterministic
    templates = sorted(f for f in pt1_dir.glob("*.in") if "adsorbed" not in f.name.lower())
    return templates


def get_frozen_positions_and_cell(cif_path, num_adsorbate_atoms=NUM_ADSORBATE_ATOMS,
                                  z_cutoff_fraction=Z_CUTOFF_FRACTION):
    """Return atom count, present elements, adsorbate element, CELL and POSITIONS cards."""
    atoms = read(cif_path)
    total_atoms = len(atoms)
    num_slab_atoms = total_atoms - num_adsorbate_atoms

    symbols = atoms.get_chemical_symbols()
    adsorbate_element = symbols[-1]
    unique_elements = list(dict.fromkeys(symbols))  # order-preserving unique

    # Slab CIFs wrap across the periodic boundary, so remap fractional z
    scaled_z = atoms.get_scaled_positions()[:, 2]
    mapped_z = (scaled_z + 0.1) % 1.0 - 0.1

    slab_z = mapped_z[:num_slab_atoms]
    z_min, z_max = slab_z.min(), slab_z.max()
    z_threshold = z_min + z_cutoff_fraction * (z_max - z_min)

    fixed_indices = {i for i in range(num_slab_atoms) if mapped_z[i] < z_threshold}

    cell_lines = ["CELL_PARAMETERS angstrom\n"]
    for vec in atoms.get_cell():
        cell_lines.append(f"  {vec[0]:14.10f}  {vec[1]:14.10f}  {vec[2]:14.10f}\n")

    pos_lines = ["ATOMIC_POSITIONS angstrom\n"]
    for i, atom in enumerate(atoms):
        x, y, z = atom.position
        flags = "   0 0 0" if i in fixed_indices else ""
        pos_lines.append(f"{atom.symbol:<2s}  {x:14.10f}  {y:14.10f}  {z:14.10f}{flags}\n")

    return total_atoms, unique_elements, adsorbate_element, "".join(cell_lines), "".join(pos_lines)


def extract_template_starting_mag(template_content):
    """Extract starting_magnetization value from the template input file.
    
    Returns the last starting_magnetization value found, or None if not present.
    """
    matches = re.findall(
        r"starting_magnetization(?:\(\s*[\d\w]+\s*\))?\s*=\s*([\d\.\-+eE]+)",
        template_content,
        flags=re.IGNORECASE
    )
    return float(matches[-1]) if matches else None


def build_atomic_species_block(elements):
    """Generate ATOMIC_SPECIES block based on elements present in the structure."""
    lines = ["ATOMIC_SPECIES\n"]
    for elem in elements:
        if elem not in PSEUDO_DB:
            raise KeyError(f"Element '{elem}' missing in PSEUDO_DB configuration.")
        info = PSEUDO_DB[elem]
        lines.append(f"  {elem:<2s}  {info['mass']:7.3f}  {info['upf']}\n")
    return "".join(lines)


def process_template_header(template_content, total_atoms, elements, adsorbate_element):
    """Update nat, ntyp, set starting_magnetization for Co and Pt, and strip old cards."""
    ntyp = len(elements)
    template_mag = extract_template_starting_mag(template_content)

    # Strip every regenerated card (incl. prefix/outdir) so we can inject a
    # unique pair below -- shared scratch would otherwise clobber other variants.
    content = strip_old_cards(template_content)

    content = re.sub(r"nat\s*=\s*\d+", f"nat = {total_atoms}", content, flags=re.IGNORECASE)
    content = re.sub(r"ntyp\s*=\s*\d+", f"ntyp = {ntyp}", content, flags=re.IGNORECASE)

    mag_entries = []
    for i, elem in enumerate(elements, start=1):
        if elem == adsorbate_element:
            # Only add starting_magnetization for adsorbate if template had one
            if template_mag is not None:
                mag_entries.append(f"  starting_magnetization({i}) = {template_mag}")
        elif elem == "Co":
            mag_entries.append(f"  starting_magnetization({i}) = 0.1")

    mag_str = "\n" + "\n".join(mag_entries) if mag_entries else ""
    content = re.sub(r"(ntyp\s*=\s*\d+)", r"\1" + mag_str, content, flags=re.IGNORECASE)

    header = content.strip()
    return f"{header}\n\n{build_atomic_species_block(elements)}\n"


def main():
    parser = argparse.ArgumentParser(description="Generate QE inputs for a single Pt atom adsorbed on an LCO slab.")
    parser.add_argument(
        "--slab",
        choices=SLAB_CHOICES,
        default=DEFAULT_SLAB,
        help=f"Slab profile to process (default: {DEFAULT_SLAB})",
    )
    parser.add_argument(
        "--sites",
        nargs="+",
        choices=SITE_ORDER,
        help="Only generate inputs for these sites (default: all)",
    )
    args = parser.parse_args()

    if not PT1_DIR.exists():
        print(f"Error: Directory '{PT1_DIR}' does not exist.")
        return

    templates = find_templates_in(PT1_DIR)
    if not templates:
        print(f"Error: No template '.in' found in {PT1_DIR}")
        return

    selected_sites = set(args.sites) if args.sites else VALID_SITES
    output_name = get_slab_profile(args.slab)["output_name"]
    all_cifs = PT1_DIR.rglob(f"Adsorbed/{output_name}/*/*.cif")
    adsorbed_cifs = sorted(p for p in all_cifs if p.parent.name in selected_sites)

    if not adsorbed_cifs:
        print(
            f"No CIF files found inside {PT1_DIR}/Adsorbed/{output_name}/ "
            f"for sites: {', '.join(sorted(selected_sites))}"
        )
        return

    for template_path in templates:
        template_content = template_path.read_text()
        template_stem = template_path.stem  # e.g., "44.1_template" or "new_template"

        for cif_path in adsorbed_cifs:
            if cif_path.stat().st_size == 0:
                print(f"Skipping empty CIF file: {cif_path}")
                continue

            try:
                total_atoms, elements, adsorbate_element, cell_str, pos_str = get_frozen_positions_and_cell(cif_path)
                header_content = process_template_header(
                    template_content, total_atoms, elements, adsorbate_element
                )
            except Exception as e:
                print(f"Error processing {cif_path.name}: {e}")
                continue

            # Build variant-specific prefix/outdir so multiple runs on same machine don't clobber each other.
            unique_prefix = variant_prefix(SYSTEM_TAG, cif_path.parent.name, template_stem)
            prefix_outdir_lines = make_prefix_directive(unique_prefix, "./tmp")
            header_content = header_content.rstrip() + "\n" + prefix_outdir_lines.rstrip() + "\n"

            output_in_path = cif_path.parent / f"INPUT_{cif_path.stem}_{template_stem}.in"

            # Skip if this variant's OUTPUT_<prefix>.out already exists (idempotent runs).
            output_prefix = f"OUTPUT_{unique_prefix}"
            has_existing = any(
                child.is_file() and child.name.startswith(output_prefix) and child.suffix.lower() == ".out"
                for child in cif_path.parent.iterdir()
            )
            if has_existing:
                print(f"Skipping [{cif_path.parent.name}]: QE output already exists in {cif_path.parent}")
                continue

            with open(output_in_path, "w") as f:
                f.write(header_content + "\n")
                f.write(cell_str + "\n")
                f.write(pos_str + "\n")
                f.write(HUBBARD_BLOCK + "\n")
                f.write(SURFACE_KPOINTS)

            print(f"Generated [{cif_path.parent.name}]: {output_in_path.name} (using {template_path.name})")
if __name__ == "__main__":
    main()
