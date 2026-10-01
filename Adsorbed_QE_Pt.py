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
from adsorption_sites import VALID_SITES

# -------------------------------------------------------------------------
# Configuration & Path Resolution
# -------------------------------------------------------------------------
SCRIPT_PATH = Path(__file__).resolve()

PROJECT_ROOT = next(
    (p for p in [SCRIPT_PATH.parent] + list(SCRIPT_PATH.parents) if (p / "05_metals").exists()),
    SCRIPT_PATH.parent
)
PT1_DIR = PROJECT_ROOT / "05_metals" / "Pt1"

NUM_ADSORBATE_ATOMS = 1  # single Pt atom, always the last atom in the CIF
Z_CUTOFF_FRACTION = 0.5  # Freeze bottom 50% of slab height

# Reference pseudopotentials and atomic masses
PSEUDO_DB = {
    "Li": {"mass": 6.941,  "upf": "li_pbe_v1.4.uspp.F.UPF"},
    "Co": {"mass": 58.933, "upf": "co_pbe_v1.2.uspp.F.UPF"},
    "O":  {"mass": 15.999, "upf": "o_pbe_v1.2.uspp.F.UPF"},
    "Pt": {"mass": 195.084,"upf": "pt_pbe_v1.4.uspp.F.UPF"},
}

SURFACE_KPOINTS = "K_POINTS (automatic)\n4 3 1 0 0 0\n"
HUBBARD_BLOCK = "HUBBARD (ortho-atomic)\nU Co-3d 3.2\n"

CARDS_TO_STRIP = ("ATOMIC_SPECIES", "K_POINTS", "HUBBARD", "CELL_PARAMETERS", "ATOMIC_POSITIONS")


def find_template_in(pt1_dir):
    """Find the bare-atom template .in file directly in the Pt1 directory."""
    if not pt1_dir.exists():
        return None

    # Ignore any 'adsorbed' input files; sort so the choice is deterministic
    templates = sorted(f for f in pt1_dir.glob("*.in") if "adsorbed" not in f.name.lower())
    return templates[0] if templates else None


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
    """Extract starting_magnetization value from the template input file."""
    matches = re.findall(
        r"starting_magnetization(?:\(\s*[\d\w]+\s*\))?\s*=\s*([\d\.\-+eE]+)",
        template_content,
        flags=re.IGNORECASE
    )
    return float(matches[-1]) if matches else 0.1


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

    content = re.sub(
        r"starting_magnetization\(\s*[\d\w]+\s*\)\s*=\s*[\d\.eE+-]+,?",
        "",
        template_content,
        flags=re.IGNORECASE
    )

    content = re.sub(r"nat\s*=\s*\d+", f"nat = {total_atoms}", content, flags=re.IGNORECASE)
    content = re.sub(r"ntyp\s*=\s*\d+", f"ntyp = {ntyp}", content, flags=re.IGNORECASE)

    mag_entries = []
    for i, elem in enumerate(elements, start=1):
        if elem == adsorbate_element:
            mag_entries.append(f"  starting_magnetization({i}) = {template_mag}")
        elif elem == "Co":
            mag_entries.append(f"  starting_magnetization({i}) = 0.1")

    mag_str = "\n" + "\n".join(mag_entries) if mag_entries else ""
    content = re.sub(r"(ntyp\s*=\s*\d+)", r"\1" + mag_str, content, flags=re.IGNORECASE)

    cleaned_lines = []
    skip_block = False
    for line in content.splitlines():
        stripped = line.strip().upper()
        if stripped.startswith(CARDS_TO_STRIP):
            skip_block = True
            continue

        if skip_block and (stripped.startswith("&") or stripped.startswith("OCCUPATIONS")):
            skip_block = False

        if not skip_block and line.strip():
            cleaned_lines.append(line)

    header = "\n".join(cleaned_lines).strip()
    return f"{header}\n\n{build_atomic_species_block(elements)}\n"


def main():
    parser = argparse.ArgumentParser(description="Generate QE inputs for a single Pt atom adsorbed on an LCO slab.")
    parser.add_argument(
        "--slab",
        choices=SLAB_CHOICES,
        default=DEFAULT_SLAB,
        help=f"Slab profile to process (default: {DEFAULT_SLAB})",
    )
    args = parser.parse_args()

    if not PT1_DIR.exists():
        print(f"Error: Directory '{PT1_DIR}' does not exist.")
        return

    template_in = find_template_in(PT1_DIR)
    if template_in is None:
        print(f"Error: No template '.in' found in {PT1_DIR}")
        return

    output_name = get_slab_profile(args.slab)["output_name"]
    all_cifs = PT1_DIR.rglob(f"Adsorbed/{output_name}/*/*.cif")
    adsorbed_cifs = sorted(p for p in all_cifs if p.parent.name in VALID_SITES)

    if not adsorbed_cifs:
        print(
            f"No valid CIF files found inside {PT1_DIR}/Adsorbed/{output_name}/"
            "{li,co,lico,centered,oxygen,lili,coco}/"
        )
        return

    template_content = template_in.read_text()

    for cif_path in adsorbed_cifs:
        if cif_path.stat().st_size == 0:
            print(f"Skipping empty CIF file: {cif_path}")
            continue

        has_existing_output = any(
            child.is_file() and child.suffix.lower() == ".out"
            for child in cif_path.parent.iterdir()
        )
        if has_existing_output:
            print(f"Skipping [{cif_path.parent.name}]: QE output already exists in {cif_path.parent}")
            continue

        try:
            total_atoms, elements, adsorbate_element, cell_str, pos_str = get_frozen_positions_and_cell(cif_path)
            header_content = process_template_header(
                template_content, total_atoms, elements, adsorbate_element
            )
        except Exception as e:
            print(f"Error processing {cif_path.name}: {e}")
            continue

        output_in_path = cif_path.parent / f"INPUT_{cif_path.stem}.in"
        with open(output_in_path, "w") as f:
            f.write(header_content + "\n")
            f.write(cell_str + "\n")
            f.write(pos_str + "\n")
            f.write(HUBBARD_BLOCK + "\n")
            f.write(SURFACE_KPOINTS)

        print(f"Generated [{cif_path.parent.name}]: {output_in_path.name} (using {template_in.name})")


if __name__ == "__main__":
    main()