#!/usr/bin/env python3
"""Build single-metal-atom adsorbed VASP calculation decks on an LCO slab.

Places a metal atom (e.g. Ni, Pt) at specified adsorption sites above
a relaxed LCO slab CONTCAR/POSCAR, and generates complete VASP calculation
directories (POSCAR, INCAR, KPOINTS, POTCAR, POTCAR_SPECIES, run_vasp.sh).

Defaults to:
  --slab:     04_slab_vasp_new/CONTCAR
  --metal:    Ni
  --height:   2.50 Å
  --out-dir:  05_metals_vasp/{metal}1
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from ase import Atoms
from ase.io import read, write

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from lco.adsorption_sites import SITE_SPECS
from lco.ase_io import freeze_by_z

DEFAULT_SLAB_PATH = PROJECT_ROOT / "04_slab_vasp_new" / "CONTCAR"
DEFAULT_PP_PATH = Path(os.environ.get("VASP_PP_PATH", "/home/mark/vasp_paw_library/PBE64"))

# Default PAW potential mapping
DEFAULT_PAW_DATASETS: Dict[str, str] = {
    "Li": "Li_sv",
    "Co": "Co_pv",
    "O": "O",
    "Ni": "Ni_pv",
    "Pt": "Pt",
    "Cu": "Cu_pv",
}

# Default initial magnetic moments (uB) for adsorbate metals
DEFAULT_AD_MAGMOM: Dict[str, float] = {
    "Ni": 1.0,
    "Pt": 0.0,
    "Co": 1.0,
    "Cu": 0.0,
}

SURFACE_Z_TOLERANCE = 0.8
SURFACE_FRAC_Z_MAX = 0.8


def parse_selective_dynamics_from_vasp(file_path: Path) -> Optional[List[bool]]:
    """Parse selective dynamics boolean flags from a POSCAR/CONTCAR file.
    
    Returns a list of booleans where True means FIXED (F F F), False means FREE (T T T),
    matching lco.ase_io.freeze_by_z convention.
    """
    if not file_path.exists():
        return None
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [line.strip() for line in f if line.strip()]

    if len(lines) < 8:
        return None

    # Check if selective dynamics is present (usually line 7 or 8)
    sd_line_idx = -1
    for idx in (7, 8):
        if idx < len(lines) and lines[idx].lower().startswith("s"):
            sd_line_idx = idx
            break

    if sd_line_idx == -1:
        return None

    # Line after selective dynamics is coord format (Direct or Cartesian)
    start_idx = sd_line_idx + 2
    # Count of atoms from line 6
    try:
        atom_counts = [int(c) for c in lines[6].split()]
        total_atoms = sum(atom_counts)
    except ValueError:
        return None

    flags: List[bool] = []
    for line in lines[start_idx : start_idx + total_atoms]:
        parts = line.split()
        if len(parts) >= 6:
            # e.g., x y z T T T or F F F
            is_fixed = (parts[-3].upper() == "F")
            flags.append(is_fixed)
        else:
            flags.append(False)

    return flags if len(flags) == total_atoms else None


def get_top_surface_indices(slab: Atoms, element: str, z_tol: float = SURFACE_Z_TOLERANCE) -> List[int]:
    """Find top-most atoms of a chosen element in the slab surface layer."""
    valid_mask = slab.get_scaled_positions()[:, 2] < SURFACE_FRAC_Z_MAX
    matching = [i for i in np.where(valid_mask)[0] if slab.symbols[i] == element]
    if not matching:
        return []
    max_z = np.max(slab.positions[matching, 2])
    return [i for i in matching if slab.positions[i, 2] > (max_z - z_tol)]


def closest_to_center(slab: Atoms, indices: Sequence[int], center_xy: np.ndarray) -> int:
    """Return index of the atom whose xy coordinates are nearest the slab center."""
    return min(indices, key=lambda i: float(np.linalg.norm(slab.positions[i][:2] - center_xy)))


def surface_pairs(slab: Atoms, surface: Dict[str, List[int]], elements: Tuple[str, ...],
                  max_dist: Optional[float]) -> List[Tuple[int, int]]:
    """Generate pairs of surface atoms for bridge sites, optionally cutoff-filtered."""
    from itertools import combinations, product

    el_a, el_b = elements
    if el_a == el_b:
        pairs = list(combinations(surface[el_a], 2))
    else:
        pairs = list(product(surface[el_a], surface[el_b]))

    if max_dist is None:
        return pairs

    return [
        (i, j) for i, j in pairs
        if slab.get_distance(i, j, mic=True) < max_dist
    ]


def find_best_midpoint(slab: Atoms, pairs: List[Tuple[int, int]], center_xy: np.ndarray) -> Optional[np.ndarray]:
    """Find the midpoint of atom pairs closest to the slab xy center."""
    best_midpoint, min_dist = None, float("inf")
    for i, j in pairs:
        vec = slab.get_distance(i, j, vector=True, mic=True)
        midpoint_xy = (slab.positions[i] + vec / 2.0)[:2]
        dist = float(np.linalg.norm(midpoint_xy - center_xy))
        if dist < min_dist:
            min_dist, best_midpoint = dist, midpoint_xy
    return best_midpoint


def resolve_site(slab: Atoms, spec: dict, surface: Dict[str, List[int]], center_xy: np.ndarray) -> Optional[np.ndarray]:
    """Return xy coordinate for a site specification."""
    kind = spec["kind"]
    if kind == "center":
        return center_xy
    if kind == "atop":
        el = spec["elements"][0]
        if not surface.get(el):
            return None
        idx = closest_to_center(slab, surface[el], center_xy)
        return slab.positions[idx][:2]
    if kind == "bridge":
        pairs = surface_pairs(slab, surface, spec["elements"], spec.get("max_dist"))
        return find_best_midpoint(slab, pairs, center_xy)
    raise ValueError(f"Unknown site kind: {kind!r}")


def build_target_sites(slab: Atoms) -> Tuple[float, Dict[str, np.ndarray]]:
    """Resolve all target sites on the top surface of the slab."""
    valid_mask = slab.get_scaled_positions()[:, 2] < SURFACE_FRAC_Z_MAX
    surface_z = float(np.max(slab.positions[valid_mask, 2]))

    cell = slab.get_cell()
    center_xy = (cell[0][:2] + cell[1][:2]) / 2.0

    needed_elements = sorted({el for spec in SITE_SPECS.values() for el in spec["elements"]})
    surface = {el: get_top_surface_indices(slab, el) for el in needed_elements}

    sites: Dict[str, np.ndarray] = {}
    for name, spec in SITE_SPECS.items():
        xy = resolve_site(slab, spec, surface, center_xy)
        if xy is not None:
            sites[name] = xy
        else:
            print(f"Warning: could not resolve site {name}")

    return surface_z, sites


def generate_poscar_content(adsorbed: Atoms, selective_dynamics: List[bool], comment: str) -> str:
    """Format VASP POSCAR with Selective Dynamics."""
    cell = adsorbed.get_cell()
    symbols = list(adsorbed.get_chemical_symbols())

    species: List[str] = []
    counts: List[int] = []
    current_species = None
    for s in symbols:
        if s != current_species:
            species.append(s)
            counts.append(1)
            current_species = s
        else:
            counts[-1] += 1

    lines = [
        comment,
        "   1.00000000000000",
    ]
    for vec in cell:
        lines.append(f"     {vec[0]:21.16f}    {vec[1]:21.16f}    {vec[2]:21.16f}")
    lines.append("   " + "   ".join(species))
    lines.append("    " + "    ".join(str(c) for c in counts))
    lines.append("Selective dynamics")
    lines.append("Direct")

    scaled_positions = adsorbed.get_scaled_positions()
    for (x, y, z), is_fixed in zip(scaled_positions, selective_dynamics):
        flag = "F" if is_fixed else "T"
        lines.append(f"     {x:21.16f}    {y:21.16f}    {z:21.16f}   {flag}   {flag}   {flag}")

    return "\n".join(lines) + "\n"


def generate_incar_content(species_order: List[str], species_counts: Dict[str, int],
                           metal: str, site_name: str,
                           magmom_metal: float = 1.0) -> str:
    """Generate VASP INCAR adhering to converged slab settings with U on Co only."""
    magmom_parts: List[str] = []
    ldaul_parts: List[str] = []
    ldauu_parts: List[str] = []
    ldauj_parts: List[str] = []

    for sp in species_order:
        n = species_counts[sp]
        if sp in ("Li", "O"):
            magmom_parts.append(f"{n}*0.0")
            ldaul_parts.append("-1")
            ldauu_parts.append("0")
            ldauj_parts.append("0")
        elif sp == "Co":
            magmom_parts.append(f"{n}*1.0")
            ldaul_parts.append("2")
            ldauu_parts.append("3.2")
            ldauj_parts.append("0")
        else:
            # Other metal adatoms (e.g. Ni, Pt, Cu)
            magmom_parts.append(f"{n}*{magmom_metal:.1f}")
            ldaul_parts.append("-1")
            ldauu_parts.append("0")
            ldauj_parts.append("0")

    magmom_str = " ".join(magmom_parts)
    ldaul_str = " ".join(ldaul_parts)
    ldauu_str = " ".join(ldauu_parts)
    ldauj_str = " ".join(ldauj_parts)

    incar_lines = [
        f"SYSTEM = LiCoO2_104_slab_3x2_{metal}1_{site_name}_PBE_U3.2_D3",
        "GGA = PE",
        "ENCUT = 550",
        "PREC = Accurate",
        "EDIFF = 1E-6",
        "NELM = 500",
        "ALGO = Normal",
        "LREAL = .FALSE.",
        "LASPH = .TRUE.",
        "ISMEAR = 0",
        "SIGMA = 0.13606",
        "ISPIN = 2",
        f"MAGMOM = {magmom_str}",
        "LDAU = .TRUE.",
        "LDAUTYPE = 2",
        f"LDAUL = {ldaul_str}",
        f"LDAUU = {ldauu_str}",
        f"LDAUJ = {ldauj_str}",
        "LMAXMIX = 4",
        "IVDW = 11",
        "IBRION = 2",
        "ISIF = 2",
        "NSW = 100",
        "EDIFFG = -0.0103",
        "PSTRESS = 0",
        "LDIPOL = .TRUE.",
        "IDIPOL = 3",
        "DIPOL = 0.5 0.5 0.0935",
        "KPAR = 2",
        "NCORE = 4",
        "ISTART = 0",
        "ICHARG = 2",
        "POTIM = 0.1",
        "LORBIT = 11",
    ]
    return "\n".join(incar_lines) + "\n"


def generate_kpoints_content() -> str:
    """Generate Monkhorst-Pack Gamma-centered 4x3x1 KPOINTS."""
    return (
        "Automatic Gamma mesh for LiCoO2 (104) slab + adatom\n"
        "0\n"
        "Gamma\n"
        "4 3 1\n"
        "0 0 0\n"
    )


def generate_run_vasp_content(job_name: str) -> str:
    """Generate SLURM batch submission script matching the converged slab."""
    return f"""#!/usr/bin/env bash
# Submit from the calculation folder with: sbatch run_vasp.sh
#SBATCH --job-name={job_name}
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --output=Output.%j.out
#SBATCH --error=Error.%j.err

set -eo pipefail

if [[ -z ${{SLURM_JOB_ID:-}} || -z ${{SLURM_SUBMIT_DIR:-}} ]]; then
    echo 'Run this script through sbatch from the calculation directory.' >&2
    exit 2
fi
if [[ ${{SLURM_NTASKS:-0}} != 40 || ${{SLURM_CPUS_PER_TASK:-1}} != 1 || ${{SLURM_JOB_NUM_NODES:-${{SLURM_NNODES:-0}}}} != 1 ]]; then
    echo 'Expected one node, 40 MPI ranks, and one CPU per rank.' >&2
    exit 2
fi

cd -- "$SLURM_SUBMIT_DIR"
for input in INCAR POSCAR KPOINTS POTCAR; do
    if [[ ! -f $input || ! -r $input || ! -s $input ]]; then
        printf 'Missing, empty, or unreadable input: %s\\n' "$input" >&2
        exit 2
    fi
done

oneapi_setvars=${{ONEAPI_SETVARS:-/opt/intel/oneapi/setvars.sh}}
if [[ ! -r $oneapi_setvars ]]; then
    printf 'Cannot read oneAPI environment script: %s\\n' "$oneapi_setvars" >&2
    exit 2
fi

set +e
source "$oneapi_setvars"
setup_status=$?
set -e
if (( setup_status != 0 )); then
    printf 'oneAPI setup failed with status %s\\n' "$setup_status" >&2
    exit "$setup_status"
fi
set -u

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export I_MPI_HYDRA_BOOTSTRAP=slurm

vasp_exe=${{VASP_EXE:-/home/Software/VASP/vasp.6.5.1/bin/vasp_std}}
command -v mpiexec >/dev/null || {{ echo 'Intel MPI mpiexec not found.' >&2; exit 2; }}
command -v "$vasp_exe" >/dev/null || {{ echo 'VASP executable not found.' >&2; exit 2; }}

printf 'Job: %s\\nNodes: %s\\nMPI ranks: %s\\n' "$SLURM_JOB_ID" "${{SLURM_JOB_NODELIST:-unknown}}" "$SLURM_NTASKS"
printf 'MPI launcher: %s\\nVASP executable: %s\\n' "$(command -v mpiexec)" "$(command -v "$vasp_exe")"
printf 'Start: %s\\n' "$(date -Is)"

status=0
mpiexec -np "$SLURM_NTASKS" "$vasp_exe" || status=$?
printf 'VASP exit code: %s\\nFinished: %s\\n' "$status" "$(date -Is)"
exit "$status"
"""


def assemble_potcar(potcar_species_list: List[str], out_potcar_path: Path, pp_path: Path) -> bool:
    """Concatenate individual PAW POTCAR files into a combined POTCAR."""
    if not pp_path.exists():
        print(f"Warning: PAW path '{pp_path}' not found. Cannot assemble POTCAR.")
        return False

    chunks: List[str] = []
    for sp in potcar_species_list:
        candidate = pp_path / sp / "POTCAR"
        if not candidate.exists():
            candidate = pp_path / "PBE" / sp / "POTCAR"
        if not candidate.exists():
            print(f"Warning: POTCAR for '{sp}' not found under '{pp_path}'. Cannot assemble POTCAR.")
            return False
        text = candidate.read_text(encoding="utf-8", errors="ignore")
        # Ensure each dataset ends with exactly one newline and no blank lines between datasets
        text = text.rstrip("\n") + "\n"
        chunks.append(text)

    out_potcar_path.write_text("".join(chunks), encoding="utf-8")
    return True


def create_site_calculation(site_name: str, target_xy: np.ndarray, surface_z: float,
                            slab: Atoms, slab_selective_dynamics: List[bool],
                            metal: str, height: float, out_dir: Path, pp_path: Path,
                            paw_mapping: Dict[str, str], overwrite: bool = False) -> Path:
    """Create a single site calculation folder and write all inputs."""
    site_folder = out_dir / site_name.capitalize()
    site_folder.mkdir(parents=True, exist_ok=True)

    # 1. Adsorbed structure
    ad_atom = Atoms(metal, positions=[[target_xy[0], target_xy[1], surface_z + height]])
    raw_atoms = slab.copy()
    raw_atoms.extend(ad_atom)
    raw_sd = list(slab_selective_dynamics) + [False]  # adatom is free (T T T)

    # Group atoms and selective_dynamics flags by species preserving slab species order
    species_order: List[str] = []
    for s in slab.get_chemical_symbols():
        if s not in species_order:
            species_order.append(s)
    if metal not in species_order:
        species_order.append(metal)

    ordered_indices: List[int] = []
    all_symbols = raw_atoms.get_chemical_symbols()
    for spec in species_order:
        for idx, sym in enumerate(all_symbols):
            if sym == spec:
                ordered_indices.append(idx)

    adsorbed = raw_atoms[ordered_indices]
    selective_dynamics = [raw_sd[i] for i in ordered_indices]

    # Species counts
    species_symbols = list(adsorbed.get_chemical_symbols())
    counts: Dict[str, int] = {}
    for s in species_symbols:
        counts[s] = counts.get(s, 0) + 1

    # 2. POSCAR
    poscar_content = generate_poscar_content(
        adsorbed, selective_dynamics,
        comment=f"LiCoO2 (104) slab + {metal}1 @ {site_name} (height={height:.2f} A)"
    )
    poscar_path = site_folder / "POSCAR"
    poscar_path.write_text(poscar_content, encoding="utf-8")

    # Also save CIF for visual inspection and archiving
    cif_path = site_folder / f"adsorbed_{metal}1_{site_name}.cif"
    write(cif_path, adsorbed)

    # 3. INCAR
    ad_magmom = DEFAULT_AD_MAGMOM.get(metal, 1.0)
    incar_content = generate_incar_content(
        species_order=species_order,
        species_counts=counts,
        metal=metal,
        site_name=site_name,
        magmom_metal=ad_magmom,
    )
    (site_folder / "INCAR").write_text(incar_content, encoding="utf-8")

    # 4. KPOINTS
    kpoints_content = generate_kpoints_content()
    (site_folder / "KPOINTS").write_text(kpoints_content, encoding="utf-8")

    # 5. POTCAR_SPECIES and POTCAR
    potcar_species_list = [paw_mapping.get(s, s) for s in species_order]
    (site_folder / "POTCAR_SPECIES").write_text(" ".join(potcar_species_list) + "\n", encoding="utf-8")

    potcar_success = assemble_potcar(potcar_species_list, site_folder / "POTCAR", pp_path)

    # 6. run_vasp.sh
    run_vasp_content = generate_run_vasp_content(f"{metal}1_{site_name}")
    run_script = site_folder / "run_vasp.sh"
    run_script.write_text(run_vasp_content, encoding="utf-8")
    run_script.chmod(0o755)

    print(f"Created VASP deck for [{site_name}] in: {site_folder}")
    if potcar_success:
        print(f"  -> Assembled POTCAR ({' '.join(potcar_species_list)})")
    else:
        print(f"  -> Note: Remember to run setup_potcar.py to generate POTCAR")

    return site_folder


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build single-metal-atom adsorbed VASP decks on a relaxed LCO slab."
    )
    parser.add_argument("--slab", type=Path, default=DEFAULT_SLAB_PATH,
                        help=f"Path to relaxed slab CONTCAR/POSCAR (default: {DEFAULT_SLAB_PATH})")
    parser.add_argument("--metal", type=str, default="Ni",
                        help="Adsorbate metal element symbol (default: Ni)")
    parser.add_argument("--sites", nargs="+", default=["all"],
                        help="Sites to generate (e.g. oxygen co li lico ... or 'all')")
    parser.add_argument("--height", type=float, default=2.50,
                        help="Initial vertical adsorption distance in Angstrom (default: 2.50)")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory base (default: 05_metals_vasp/<metal>1)")
    parser.add_argument("--pp-path", type=Path, default=DEFAULT_PP_PATH,
                        help=f"Path to VASP PAW datasets (default: {DEFAULT_PP_PATH})")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing files in output folders")

    args = parser.parse_args()

    slab_path = args.slab.resolve()
    if not slab_path.exists():
        sys.exit(f"Error: Slab file does not exist: {slab_path}")

    slab = read(slab_path)
    print(f"Loaded slab from {slab_path}: {len(slab)} atoms")

    # Read selective dynamics from slab file or compute via freeze_by_z
    slab_sd = parse_selective_dynamics_from_vasp(slab_path)
    if slab_sd is not None:
        print(f"Extracted selective dynamics from slab file: {sum(slab_sd)} fixed, {len(slab_sd) - sum(slab_sd)} free")
    else:
        slab_sd = freeze_by_z(slab, z_cutoff_fraction=0.5)
        print(f"Computed selective dynamics (bottom 50% fixed): {sum(slab_sd)} fixed, {len(slab_sd) - sum(slab_sd)} free")

    # Build target adsorption sites
    surface_z, all_sites = build_target_sites(slab)
    print(f"Surface Z maximum: {surface_z:.4f} A")

    # Filter sites
    requested_sites = [s.lower() for s in args.sites]
    if "all" in requested_sites:
        target_sites = all_sites
    else:
        target_sites = {name: xy for name, xy in all_sites.items() if name.lower() in requested_sites}
        missing = [s for s in requested_sites if s not in all_sites]
        if missing:
            print(f"Warning: site(s) {missing} not recognized. Available: {list(all_sites.keys())}")

    if not target_sites:
        sys.exit("Error: No valid sites to generate.")

    # Base output folder: 05_metals_vasp/<metal>1
    metal = args.metal.capitalize()
    if args.out_dir is not None:
        base_out_dir = args.out_dir
    else:
        base_out_dir = PROJECT_ROOT / "05_metals_vasp" / f"{metal}1"

    print(f"\nGenerating {len(target_sites)} site(s) for {metal}1 in {base_out_dir}:")
    for site_name, xy in target_sites.items():
        create_site_calculation(
            site_name=site_name,
            target_xy=xy,
            surface_z=surface_z,
            slab=slab,
            slab_selective_dynamics=slab_sd,
            metal=metal,
            height=args.height,
            out_dir=base_out_dir,
            pp_path=args.pp_path,
            paw_mapping=DEFAULT_PAW_DATASETS,
            overwrite=args.overwrite,
        )

    print("\nAll VASP input decks created successfully!")


if __name__ == "__main__":
    main()

