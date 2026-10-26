"""Repository path resolution.

Every script in this project used to resolve the repository root on its own,
and one of them hard-coded an absolute ``/home/mark/...`` path.  This module
is the single source of truth: :func:`project_root` walks upwards from a
starting directory until it finds the marker directory that identifies the
repository, so the code works no matter where it is invoked from.
"""
from pathlib import Path

#: Directory names that identify the repository root.
ROOT_MARKERS = ("05_metals", "04_slab")


def project_root(start=None):
    """Return the repository root.

    The search starts at ``start`` (default: this file's directory) and walks
    upwards, returning the first ancestor that contains a marker directory.
    If no marker is found the starting directory is returned, so callers always
    get a usable path.
    """
    base = Path(start).resolve() if start is not None else Path(__file__).resolve().parent
    candidates = [base, *base.parents]
    for candidate in candidates:
        if any((candidate / marker).is_dir() for marker in ROOT_MARKERS):
            return candidate
    return base


PROJECT_ROOT = project_root()
SLAB_DIR = PROJECT_ROOT / "04_slab"
METALS_DIR = PROJECT_ROOT / "05_metals"


def slab_profile_dir(profile_name):
    """Return the ``04_slab`` directory holding a named slab profile."""
    return SLAB_DIR / profile_name
