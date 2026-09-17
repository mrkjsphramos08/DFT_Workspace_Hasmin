"""Archive the legacy 44.1/441 cluster runs out of ``05_metals``.

The cluster directories under ``05_metals`` mix two generations of Quantum
ESPRESSO calculations:

* the original ``ecutwfc = 44.1`` / ``ecutrho = 441`` runs, whose files carry
  no cutoff tag in their name, and
* the re-run ``ecutwfc = 40.4`` / ``ecutrho = 500`` runs, tagged with
  ``_ecut40.4_500`` in the file name.

Only the tagged files are current. This script moves the untagged generation
to a sibling directory (``05_metals_legacy`` by default) so what remains in
``05_metals`` is unambiguous. The single most important consequence: the
legacy ``INPUT_<metal>_<shape>.in`` template disappears, so
``Adsorbed_QE_lico.py`` -- which takes the first non-adsorbed ``*.in`` it
globs as its template -- can no longer pick up the obsolete 44.1/441 file.

Safety notes
------------
* ``tmp/`` scratch sitting directly in ``<cluster>/<shape>/`` belongs to the
  *current* tagged runs and is never touched. (Old and new inputs share the
  same ``prefix``, so the tagged runs reused those directories in place.)
* Stale scratch left behind by the old adsorbed runs lives under
  ``Adsorbed/``; it is only removed when ``--tmp-mode delete`` is requested.
  Every candidate is checked against a staleness guard before removal.
* Structure files (``*_latest.cif``), ``run.sh`` and the ``Adsorbed/<site>``
  CIFs stay where they are.

Usage::

    python archive_legacy_cluster_runs.py                    # dry run
    python archive_legacy_cluster_runs.py --tmp-mode delete   # dry run + scratch plan
    python archive_legacy_cluster_runs.py --apply --tmp-mode delete
"""

import argparse
import os
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
METALS_DIR = PROJECT_ROOT / "05_metals"
DEFAULT_DEST = PROJECT_ROOT / "05_metals_legacy"

# Marks the current 40.4/500 generation.
NEW_TAG = "_ecut"

SHAPE_NAMES = ("Tetra", "Rhom")

# Placeholder directories that never received a calculation.
STUB_CLUSTERS = ("PtCo", "PtCu", "PtNi")

# Scratch directories that belong to the current tagged cluster runs.
# Stored as glob patterns so they can be passed straight to Path.glob().
PROTECTED_SCRATCH = tuple(
    f"*/{shape_name}/tmp" for shape_name in SHAPE_NAMES
)


def path_size(path):
    """Return the size in bytes of a file, or of a directory tree."""
    if path.is_file():
        return path.stat().st_size
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            file_path = Path(root) / name
            try:
                total += file_path.stat().st_size
            except OSError:
                continue
    return total


def human_size(num_bytes):
    """Format a byte count using binary units."""
    value = float(num_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{int(value)} B" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{value:.1f} TiB"


def newest_tagged_output(cluster_dir):
    """Return the newest tagged output in one cluster family, if any."""
    tagged = [
        path
        for path in cluster_dir.glob(f"*/OUTPUT_*{NEW_TAG}*.out")
        if path.is_file()
    ]
    if not tagged:
        return None
    return max(tagged, key=lambda path: path.stat().st_mtime)


def newest_mtime(path):
    """Return the newest modification time anywhere inside ``path``."""
    newest = path.stat().st_mtime
    if path.is_dir():
        for root, _dirs, files in os.walk(path):
            for name in files:
                try:
                    newest = max(newest, (Path(root) / name).stat().st_mtime)
                except OSError:
                    continue
    return newest


def collect_moves():
    """Return ``(source, reason)`` pairs for every legacy item to relocate."""
    items = []

    # 1. Legacy cluster-level input/output pairs, e.g. INPUT_Pt4_tetra.in.
    for shape_name in SHAPE_NAMES:
        for cluster_dir in sorted(METALS_DIR.glob(f"*/{shape_name}")):
            for pattern in ("INPUT_*.in", "OUTPUT_*.out"):
                for path in sorted(cluster_dir.glob(pattern)):
                    if NEW_TAG in path.stem:
                        continue
                    items.append((path, f"legacy {shape_name} input/output"))

    # 2. Legacy adsorbed-run inputs, outputs and SLURM logs.
    adsorbed_patterns = {
        "*/Tetra/Adsorbed/*/INPUT_adsorbed_*.in": "legacy adsorbed input",
        "*/Tetra/Adsorbed/*/OUTPUT_adsorbed_*.out": "legacy adsorbed output",
        "*/Tetra/Adsorbed/*/slurm-*.out": "legacy adsorbed SLURM log",
        "*/Tetra/Adsorbed/INPUT_adsorbed_*.in": "legacy orphaned adsorbed input",
        "*/Tetra/Adsorbed/OUTPUT_adsorbed_*.out": "legacy orphaned adsorbed output",
    }
    for pattern, reason in adsorbed_patterns.items():
        for path in sorted(METALS_DIR.glob(pattern)):
            items.append((path, reason))

    # 3. Empty placeholder cluster directories.
    for name in STUB_CLUSTERS:
        stub_dir = METALS_DIR / name
        if stub_dir.is_dir():
            items.append((stub_dir, "empty cluster stub"))

    return items


def collect_scratch(moves):
    """Return the stale scratch directories that ``--tmp-mode delete`` removes."""
    moved = {source for source, _reason in moves}

    candidates = set()
    for pattern in ("*/Tetra/Adsorbed/tmp", "*/Tetra/Adsorbed/*/tmp"):
        candidates.update(METALS_DIR.glob(pattern))
    for name in STUB_CLUSTERS:
        candidates.add(METALS_DIR / name / "tmp")

    protected = set()
    for pattern in PROTECTED_SCRATCH:
        protected.update(METALS_DIR.glob(pattern))

    stale = []
    skipped = [
        (path, "current tagged-run scratch (protected)")
        for path in sorted(protected)
        if path.is_dir()
    ]
    for scratch in sorted(candidates):
        if not scratch.is_dir():
            continue
        if scratch in protected:
            continue
        if any(scratch == source or scratch.is_relative_to(source) for source in moved):
            continue

        # Staleness guard: refuse to delete scratch that is newer than the
        # newest tagged output of the same cluster family.
        cluster_dir = METALS_DIR / scratch.relative_to(METALS_DIR).parts[0]
        reference = newest_tagged_output(cluster_dir)
        if reference is not None and newest_mtime(scratch) > reference.stat().st_mtime:
            skipped.append((scratch, f"newer than {reference.name} (protected)"))
            continue

        stale.append(scratch)

    return stale, skipped


def relocate(source, dest_root, apply_changes):
    """Move one item into the mirror location under ``dest_root``."""
    relative = source.relative_to(METALS_DIR)
    target = dest_root / relative

    if target.exists():
        print(f"  SKIP  {relative} -> target already exists at {target}")
        return None

    if apply_changes:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(target))
        print(f"  moved {relative}")
    else:
        print(f"  would move {relative}")
    return target


def delete_scratch(scratch, apply_changes):
    """Remove one stale scratch directory."""
    relative = scratch.relative_to(METALS_DIR)
    size = path_size(scratch)
    if apply_changes:
        shutil.rmtree(scratch)
        print(f"  deleted {relative} ({human_size(size)})")
    else:
        print(f"  would delete {relative} ({human_size(size)})")
    return size


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Move the legacy 44.1/441 cluster inputs and outputs out of 05_metals."
        )
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=DEFAULT_DEST,
        help=f"Mirror destination directory (default: {DEFAULT_DEST.name}).",
    )
    parser.add_argument(
        "--tmp-mode",
        choices=("leave", "delete"),
        default="leave",
        help=(
            "How to treat stale tmp/ scratch from the old adsorbed runs "
            "(default: leave). Current tagged-run scratch is never touched."
        ),
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Perform the moves. Without this flag, only print the planned changes.",
    )
    args = parser.parse_args()

    if not METALS_DIR.exists():
        raise FileNotFoundError(f"Metals directory missing: {METALS_DIR}")

    dest_root = args.dest.resolve()
    if dest_root == METALS_DIR or dest_root.is_relative_to(METALS_DIR):
        raise ValueError(
            f"Destination {dest_root} must live outside {METALS_DIR} so that "
            "recursive OUTPUT_*.out / *_latest.cif globs stay clean."
        )

    moves = collect_moves()
    scratch, protected_scratch = collect_scratch(moves)

    action = "Moving" if args.apply else "Would move"
    print(f"{action} {len(moves)} legacy item(s) to {dest_root}")

    moved_bytes = 0
    current_reason = None
    for source, reason in moves:
        if reason != current_reason:
            current_reason = reason
            print(f"\n[{reason}]")
        moved_bytes += path_size(source)
        relocate(source, dest_root, args.apply)

    print(f"\nTotal legacy data: {human_size(moved_bytes)}")

    if args.tmp_mode == "delete":
        deleted_bytes = sum(path_size(path) for path in scratch)
        action = "Deleting" if args.apply else "Would delete"
        print(
            f"\n{action} {len(scratch)} stale scratch dir(s): "
            f"{human_size(deleted_bytes)}"
        )
        for path in scratch:
            delete_scratch(path, args.apply)
    else:
        print(
            "\nStale scratch left in place. Re-run with --tmp-mode delete to "
            "remove it."
        )

    if protected_scratch:
        print("\nKept (belongs to the current 40.4/500 runs):")
        for path, reason in protected_scratch:
            print(f"  {path.relative_to(METALS_DIR)}: {reason}")

    if not args.apply:
        print("\nDry run only. Re-run with --apply to perform these changes.")


if __name__ == "__main__":
    main()