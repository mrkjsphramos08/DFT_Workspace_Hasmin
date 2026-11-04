"""Move legacy unlabelled adsorption folders into the ``4_layers`` namespace.

Dry-run by default; pass ``--apply`` to perform the moves.
"""
import argparse
import shutil
from pathlib import Path

from adsorption_config import PROJECT_ROOT
from adsorption_sites import VALID_SITES


METALS_DIR = PROJECT_ROOT / "05_metals"


def migrate_site_folder(site_dir, apply_changes):
    target_dir = site_dir.parent / "4_layers" / site_dir.name
    if target_dir.exists():
        print(f"Skipping {site_dir}: target already exists at {target_dir}")
        return

    action = "Moving" if apply_changes else "Would move"
    print(f"{action}: {site_dir} -> {target_dir}")
    if apply_changes:
        target_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(site_dir), str(target_dir))


def main():
    parser = argparse.ArgumentParser(
        description="Move legacy unlabelled adsorption folders into the 4_layers namespace."
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Perform the moves. Without this flag, only print the planned changes.",
    )
    args = parser.parse_args()

    if not METALS_DIR.exists():
        raise FileNotFoundError(f"Metals directory missing: {METALS_DIR}")

    legacy_dirs = [
        path
        for tetra_dir in METALS_DIR.glob("*/Tetra")
        for path in tetra_dir.glob("Adsorbed/*")
        if path.is_dir() and path.name in VALID_SITES
    ]

    if not legacy_dirs:
        print("No legacy unlabelled adsorption folders found.")
        return

    for site_dir in sorted(legacy_dirs):
        migrate_site_folder(site_dir, args.apply)


if __name__ == "__main__":
    main()