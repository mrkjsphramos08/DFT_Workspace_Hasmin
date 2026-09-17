"""Generate pw.x inputs for the adsorbed cluster systems.

Uses the bare-cluster input in the parent ``Tetra`` directory as a template,
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
METALS_DIR = PROJECT_ROOT / "05_metals"

NUM_CLUSTER_ATOMS = 4
Z_CUTOFF_FRACTION = 0.5  # Freeze bottom 50% of slab height

# Reference pseudopotentials and atomic masses
PSEUDO_DB = {
    "Li": {"mass": 6.941,  "upf": "li_pbe_v1.4.uspp.F.UPF"},
    "Co": {"mass": 58.933, "upf": "co_pbe_v1.2.uspp.F.UPF"},
    "O":  {"mass": 15.999, "upf": "o_pbe_v1.2.uspp.F.UPF"},
    "Pt": {"mass": 195.084,"upf": "pt_pbe_v1.4.uspp.F.UPF"},
    "Ni": {"mass": 58.693, "upf": "ni_pbe_v1.4.uspp.F.UPF"},
    "Cu": {"mass": 63.546, "upf": "cu_pbe_v1.2.uspp.F.UPF"},
}

SURFACE_KPOINTS = "K_POINTS (automatic)\n4 3 1 0 0 0\n"
HUBBARD_BLOCK = "HUBBARD (ortho-atomic)\nU Co-3d 3.2\n"

def find_template_in(tetra_dir):
    """Finds the template .in file directly in the Tetra directory (e.g., INPUT_Pt4_tetra.in)."""
    if not tetra_dir or not tetra_dir.exists():
        return None

    # Search directly in the Tetra directory, ignoring any 'adsorbed' input files
    templates = [
        f for f in tetra_dir.glob("*.in")
        if "adsorbed" not in f.name.lower()
    ]
    if templates:
        return templates[0]

    return None

def get_frozen_positions_and_cell(cif_path, num_cluster_atoms=4, z_cutoff_fraction=0.5):
    """Returns cell parameters, fixed positions, total atom count, present elements, and cluster element."""
    atoms = read(cif_path)
    total_atoms = len(atoms)
    num_slab_atoms = total_atoms - num_cluster_atoms

    cluster_metal = atoms.get_chemical_symbols()[-1]

    raw_symbols = list(atoms.get_chemical_symbols())
    unique_elements = []
    for s in raw_symbols:
        if s not in unique_elements:
            unique_elements.append(s)

    scaled_z = atoms.get_scaled_positions()[:, 2]
    mapped_z = (scaled_z + 0.1) % 1.0 - 0.1

    slab_mapped_z = mapped_z[:num_slab_atoms]
    z_min, z_max = slab_mapped_z.min(), slab_mapped_z.max()
    z_threshold = z_min + z_cutoff_fraction * (z_max - z_min)

    fixed_indices = set(i for i in range(num_slab_atoms) if mapped_z[i] < z_threshold)

    cell_lines = ["CELL_PARAMETERS angstrom\n"]
    for vec in atoms.get_cell():
        cell_lines.append(f"  {vec[0]:14.10f}  {vec[1]:14.10f}  {vec[2]:14.10f}\n")

    pos_lines = ["ATOMIC_POSITIONS angstrom\n"]
    for i, atom in enumerate(atoms):
        x, y, z = atom.position
        if i in fixed_indices:
            pos_lines.append(f"{atom.symbol:<2s}  {x:14.10f}  {y:14.10f}  {z:14.10f}   0 0 0\n")
        else:
            pos_lines.append(f"{atom.symbol:<2s}  {x:14.10f}  {y:14.10f}  {z:14.10f}\n")

    return total_atoms, unique_elements, cluster_metal, "".join(cell_lines), "".join(pos_lines)

def extract_template_starting_mag(template_content):
    """Extracts starting_magnetization value from the template input file."""
    matches = re.findall(
        r"starting_magnetization(?:\(\s*[\d\w]+\s*\))?\s*=\s*([\d\.\-+eE]+)",
        template_content,
        flags=re.IGNORECASE
    )
    if matches:
        return float(matches[-1])
    return 0.1

def build_atomic_species_block(elements):
    """Generates ATOMIC_SPECIES block based on elements present in the structure."""
    lines = ["ATOMIC_SPECIES\n"]
    for elem in elements:
        if elem not in PSEUDO_DB:
            raise KeyError(f"Element '{elem}' missing in PSEUDO_DB configuration.")
        info = PSEUDO_DB[elem]
        lines.append(f"  {elem:<2s}  {info['mass']:7.3f}  {info['upf']}\n")
    return "".join(lines)

def process_template_header(template_content, total_atoms, elements, cluster_metal):
    """Updates nat, ntyp, assigns starting_magnetization for Co and metal, and cleans cards."""
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
    for i, elem in enumerate(elements):
        idx = i + 1
        if elem == "Co" and elem != cluster_metal:
            mag_entries.append(f"  starting_magnetization({idx}) = 0.1")
        elif elem == cluster_metal:
            mag_entries.append(f"  starting_magnetization({idx}) = {template_mag}")

    mag_str = "\n" + "\n".join(mag_entries) if mag_entries else ""
    content = re.sub(r"(ntyp\s*=\s*\d+)", r"\1" + mag_str, content, flags=re.IGNORECASE)

    lines = content.splitlines()
    cleaned_lines = []
    skip_block = False

    cards_to_strip = [
        "ATOMIC_SPECIES", 
        "K_POINTS", 
        "HUBBARD", 
        "CELL_PARAMETERS", 
        "ATOMIC_POSITIONS"
    ]

    for line in lines:
        stripped = line.strip().upper()
        if any(stripped.startswith(card) for card in cards_to_strip):
            skip_block = True
            continue
        
        if skip_block and (stripped.startswith("&") or stripped.startswith("OCCUPATIONS")):
            skip_block = False

        if not skip_block and line.strip():
            cleaned_lines.append(line)

    header = "\n".join(cleaned_lines).strip()
    species_str = build_atomic_species_block(elements)
    
    return f"{header}\n\n{species_str}\n"

def main():
    parser = argparse.ArgumentParser(description="Generate QE inputs for adsorbed tetrahedral clusters.")
    parser.add_argument(
        "--slab",
        choices=SLAB_CHOICES,
        default=DEFAULT_SLAB,
        help=f"Slab profile to process (default: {DEFAULT_SLAB})",
    )
    args = parser.parse_args()

    if not METALS_DIR.exists():
        print(f"Error: Directory '{METALS_DIR}' does not exist relative to project root.")
        return

    output_name = get_slab_profile(args.slab)["output_name"]
    all_cifs = list(METALS_DIR.rglob(f"Tetra/Adsorbed/{output_name}/*/*.cif"))
    adsorbed_cifs = [p for p in all_cifs if p.parent.name in VALID_SITES]

    if not adsorbed_cifs:
        print(
            f"No valid CIF files found inside {METALS_DIR}/*/Tetra/Adsorbed/{output_name}/"
            "{{li,co,lico,centered,oxygen,lili,coco}}/"
        )
        return

    for cif_path in adsorbed_cifs:
        if not cif_path.exists() or cif_path.stat().st_size == 0:
            print(f"Skipping empty or missing CIF file: {cif_path}")
            continue

        tetra_dir = next((p for p in cif_path.parents if p.name == "Tetra"), None)
        if tetra_dir is None:
            print(f"Skipping {cif_path.name}: Could not locate parent 'Tetra' directory.")
            continue

        template_in = find_template_in(tetra_dir)
        if template_in is None:
            print(f"Skipping {cif_path.name}: No template '.in' found in {tetra_dir}")
            continue

        output_in_path = cif_path.parent / f"INPUT_{cif_path.stem}.in"
        has_existing_output = any(
            child.is_file() and child.suffix.lower() == ".out"
            for child in cif_path.parent.iterdir()
        )
        if has_existing_output:
            print(f"Skipping [{cif_path.parent.name}]: QE output already exists in {cif_path.parent}")
            continue

        try:
            total_atoms, elements, cluster_metal, cell_str, pos_str = get_frozen_positions_and_cell(
                cif_path, NUM_CLUSTER_ATOMS, Z_CUTOFF_FRACTION
            )
        except Exception as e:
            print(f"Error processing {cif_path.name}: {e}")
            continue

        with open(template_in, "r") as f:
            template_content = f.read()

        header_content = process_template_header(
            template_content, total_atoms, elements, cluster_metal
        )

        with open(output_in_path, "w") as f:
            f.write(header_content + "\n")
            f.write(cell_str + "\n")
            f.write(pos_str + "\n")
            f.write(HUBBARD_BLOCK + "\n")
            f.write(SURFACE_KPOINTS)

        print(f"Generated [{cif_path.parent.name}]: {output_in_path.name} (using {template_in.name})")

if __name__ == "__main__":
    main()