"""Regenerate cluster inputs at new cutoffs from a relaxed QE output.

Reads the final coordinates from an ``OUTPUT_*.out``, substitutes them into the
input template, rewrites ``ecutwfc``/``ecutrho``/``calculation`` and writes a
cutoff-tagged copy such as ``INPUT_Pt4_tetra_ecut40.4_500.in``.
"""
import argparse
import re
from pathlib import Path

from ase.io import read

PROJECT_ROOT = Path(__file__).resolve().parent
METALS_DIR = PROJECT_ROOT / "05_metals"

# Trailing cutoff tag written by an earlier run of this script, e.g. the
# "_ecut40.4_500" suffix in "INPUT_Pt4_tetra_ecut40.4_500". Existing inputs
# must be accepted as templates so the script keeps working now that the
# legacy untagged generation lives in 05_metals_legacy.
CUTOFF_TAG_RE = re.compile(r"_ecut[0-9.]+_[0-9.]+$")


def base_stem(input_path):
    """Strip any pre-existing cutoff tag so re-runs stay single-tagged."""
    return CUTOFF_TAG_RE.sub("", input_path.stem)


def resolve_output(cluster_dir, input_path, override=None):
    """Locate the QE output that holds the coordinates to reuse.

    The legacy untagged ``OUTPUT_<metal>_<shape>.out`` is preferred while it is
    still present. Once it has been archived to ``05_metals_legacy`` the newest
    cutoff-tagged output for the same cluster is used instead, so the rerun
    workflow survives the archival of the 44.1/441 generation.
    """
    if override is not None:
        return override if override.exists() else None

    legacy = cluster_dir / input_path.name.replace("INPUT_", "OUTPUT_").replace(".in", ".out")
    if legacy.exists():
        return legacy

    stem = input_path.name.replace("INPUT_", "OUTPUT_").replace(".in", "").rsplit("_ecut", 1)[0]
    candidates = sorted(
        cluster_dir.glob(f"{stem}_ecut*.out"),
        key=lambda path: path.stat().st_mtime,
    )
    return candidates[-1] if candidates else None


def format_positions(atoms):
    lines = []
    for symbol, position in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        x, y, z = position
        lines.append(f"{symbol:<2s}  {x:14.10f}  {y:14.10f}  {z:14.10f}\n")
    return lines


def update_input(input_path, output_path, ecutwfc, ecutrho, calculation, overwrite=False):
    atoms = read(output_path, format="espresso-out", index=-1)
    lines = input_path.read_text().splitlines(keepends=True)

    position_start = next(
        (index for index, line in enumerate(lines)
         if line.strip().upper().startswith("ATOMIC_POSITIONS")),
        None,
    )
    if position_start is None:
        raise ValueError(f"ATOMIC_POSITIONS block not found in {input_path}")

    position_lines = format_positions(atoms)
    replace_indices = []
    for index in range(position_start + 1, len(lines)):
        if lines[index].strip():
            replace_indices.append(index)
        if len(replace_indices) == len(position_lines):
            break
    if len(replace_indices) != len(position_lines):
        raise ValueError(f"Could not replace all atom positions in {input_path}")

    for index, position_line in zip(replace_indices, position_lines):
        lines[index] = position_line

    content = "".join(lines)
    content, wfc_count = re.subn(
        r"(\becutwfc\s*=\s*)[0-9.+-EeDd]+",
        rf"\g<1>{ecutwfc:g}",
        content,
        count=1,
        flags=re.IGNORECASE,
    )
    content, rho_count = re.subn(
        r"(\becutrho\s*=\s*)[0-9.+-EeDd]+",
        rf"\g<1>{ecutrho:g}",
        content,
        count=1,
        flags=re.IGNORECASE,
    )
    if wfc_count != 1 or rho_count != 1:
        raise ValueError(f"Could not update both cutoffs in {input_path}")

    content, calculation_count = re.subn(
        r"(\bcalculation\s*=\s*)['\"][^'\"]+['\"]",
        rf"\g<1>'{calculation}'",
        content,
        count=1,
        flags=re.IGNORECASE,
    )
    if calculation_count != 1:
        raise ValueError(f"Could not update calculation type in {input_path}")

    output_input = input_path.with_name(
        f"{base_stem(input_path)}_ecut{ecutwfc:g}_{ecutrho:g}.in"
    )
    if output_input.exists() and not overwrite:
        print(f"Skipping {output_input.relative_to(PROJECT_ROOT)}: already exists")
        return None
    output_input.write_text(content)
    print(f"Generated {output_input.relative_to(PROJECT_ROOT)}")
    return output_input


def main():
    parser = argparse.ArgumentParser(
        description="Create cluster inputs using final coordinates from QE output files."
    )
    parser.add_argument("--ecutwfc", type=float, required=True)
    parser.add_argument("--ecutrho", type=float, required=True)
    parser.add_argument(
        "--calculation",
        choices=("scf", "relax"),
        default="scf",
        help="Calculation type for the new input (default: scf).",
    )
    parser.add_argument(
        "--source-output",
        type=Path,
        help=(
            "Explicit QE output to read final coordinates from. Default: the "
            "legacy untagged OUTPUT_<metal>_<shape>.out when present, otherwise "
            "the newest OUTPUT_<metal>_<shape>_ecut*.out."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite the generated input if it already exists.",
    )
    parser.add_argument(
        "--cluster",
        help="Cluster directory name, such as Pt4. Default: process all clusters.",
    )
    parser.add_argument(
        "--shape",
        choices=("Tetra", "Rhom", "all"),
        default="all",
        help="Cluster shape to process (default: all).",
    )
    args = parser.parse_args()

    shape_names = ("Tetra", "Rhom") if args.shape == "all" else (args.shape,)
    cluster_dirs = sorted(
        path
        for shape_name in shape_names
        for path in METALS_DIR.glob(f"*/{shape_name}")
        if path.is_dir() and (args.cluster is None or path.parent.name == args.cluster)
    )
    if not cluster_dirs:
        raise FileNotFoundError("No matching cluster directories found")

    processed = 0
    for cluster_dir in cluster_dirs:
        for input_path in sorted(cluster_dir.glob("INPUT_*.in")):
            if "adsorbed" in input_path.name.lower():
                continue
            output_path = resolve_output(cluster_dir, input_path, args.source_output)
            if output_path is None:
                print(f"Skipping {input_path.name}: no QE output found to read coordinates from")
                continue
            try:
                if update_input(
                    input_path,
                    output_path,
                    args.ecutwfc,
                    args.ecutrho,
                    args.calculation,
                    overwrite=args.overwrite,
                ):
                    processed += 1
            except Exception as error:
                print(f"Error processing {input_path}: {error}")

    if processed == 0:
        print("No new inputs generated.")


if __name__ == "__main__":
    main()
