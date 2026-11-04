"""Tidy the leftover odds and ends in the ``05_metals`` adsorption folders.

After the legacy ``Adsorbed/<site>/`` layout was migrated into the ``4_layers``
namespace there are three kinds of residue left in the tree:

* orphaned ``*.cif`` files sitting directly in an ``Adsorbed`` folder instead of
  a site subdirectory (duplicates of a real site structure),
* zero-byte placeholder files, and
* ``run.sh`` scripts whose referenced ``INPUT_*.in`` no longer exists -- these
  are survivors of the archived 44.1/441 generation and some of them even
  append (``>>``) to their output, so re-running one would corrupt results.

Staleness of a ``run.sh`` is decided by its *content* (does the input file it
names exist next to it?), never by its path. That way the freshly generated
``3_layers/<site>/run.sh`` scripts can never be selected by mistake.

Items are archived rather than deleted, so the removal stays reversible.

Usage::

    python cleanup_adsorbed_folders.py            # dry run
    python cleanup_adsorbed_folders.py --apply
"""

import argparse
import re
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
METALS_DIR = PROJECT_ROOT / "05_metals"
DEFAULT_DEST = PROJECT_ROOT / "05_metals_legacy"

# Matches the input file named on the pw.x line of a run.sh.
INPUT_REF_RE = re.compile(r"pw\.x\s+-i\s+(\S+)")


def find_orphan_cifs():
    """``*.cif`` files sitting directly in an Adsorbed folder."""
    return sorted(METALS_DIR.glob("*/Tetra/Adsorbed/*.cif"))


def find_empty_files():
    """Zero-byte regular files anywhere under an Adsorbed folder."""
    return sorted(
        path
        for path in METALS_DIR.glob("*/Tetra/Adsorbed/**/*")
        if path.is_file() and path.stat().st_size == 0
    )


def find_stale_run_scripts():
    """``run.sh`` files whose referenced INPUT_*.in is missing next to them."""
    stale = []
    for run_sh in sorted(METALS_DIR.glob("*/Tetra/Adsorbed/**/run.sh")):
        match = INPUT_REF_RE.search(run_sh.read_text(errors="ignore"))
        if match is None or not (run_sh.parent / match.group(1)).exists():
            stale.append(run_sh)
    return stale


def collect_targets():
    """Return ``(path, reason)`` pairs grouped by reason, path-ordered."""
    groups = (
        ("stale run.sh (referenced INPUT_*.in is missing)", find_stale_run_scripts()),
        ("orphaned CIF outside a site directory", find_orphan_cifs()),
        ("zero-byte leftover", find_empty_files()),
    )
    seen = set()
    targets = []
    for reason, paths in groups:
        for path in paths:
            if path in seen:
                continue
            seen.add(path)
            targets.append((path, reason))
    return targets


def human_size(num_bytes):
    value = float(num_bytes)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024.0 or unit == "GiB":
            return f"{int(value)} B" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{value:.1f} GiB"


def archive(source, dest_root, apply_changes):
    """Move one leftover into the mirror location under ``dest_root``."""
    relative = source.relative_to(METALS_DIR)
    target = dest_root / relative

    if target.exists():
        print(f"  SKIP   {relative} (already archived)")
        return False

    if apply_changes:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(target))
        print(f"  moved  {relative}")
    else:
        print(f"  would  {relative}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Tidy leftover files in the 05_metals adsorption folders."
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=DEFAULT_DEST,
        help=f"Archive destination (default: {DEFAULT_DEST.name}).",
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
        raise ValueError(f"Destination {dest_root} must live outside {METALS_DIR}.")

    targets = collect_targets()
    action = "Archiving" if args.apply else "Would archive"
    print(f"{action} {len(targets)} leftover item(s) to {dest_root}\n")

    total = 0
    current_reason = None
    archived = 0
    for path, reason in targets:
        if reason != current_reason:
            current_reason = reason
            print(f"[{reason}]")
        total += path.stat().st_size
        if archive(path, dest_root, args.apply):
            archived += 1
        print()

    print(f"Total size: {human_size(total)}")
    print(f"Items archived: {archived} of {len(targets)}")

    if not args.apply:
        print("\nDry run only. Re-run with --apply to perform these changes.")

    # Report what is deliberately left in place.
    stale_set = {path for path, _reason in targets}
    kept = [
        path
        for path in sorted(METALS_DIR.glob("*/Tetra/Adsorbed/*/*/run.sh"))
        if path not in stale_set
    ]
    if kept:
        print(f"\nKept {len(kept)} active run.sh script(s), including:")
        for path in kept[:3]:
            print(f"  {path.relative_to(METALS_DIR)}")
        if len(kept) > 3:
            print(f"  ... and {len(kept) - 3} more")


if __name__ == "__main__":
    main()