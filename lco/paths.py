"""Repository path resolution.

Every script in this project used to resolve the repository root on its own,
and one of them hard-coded an absolute ``/home/mark/...`` path.  This module
is the single source of truth: :func:`project_root` walks upwards from a
starting directory until it finds the marker directory that identifies the
repository, so the code works no matter where it is invoked from.
"""
from pathlib import Path

#: Directory names that identify the repository root.
ROOT_MARKERS = (
    "05_metals_qe",
    "04_slab_qe",
    "04_slab_vasp",
    "05_metals",
    "04_slab",
)


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

# Canonical Stage directories (Quantum ESPRESSO)
BULK_QE_DIR = PROJECT_ROOT / "01_bulk_qe"
SLAB_QE_DIR = PROJECT_ROOT / "04_slab_qe"
METALS_QE_DIR = PROJECT_ROOT / "05_metals_qe"

# Canonical Stage directories (VASP)
BULK_VASP_DIR = PROJECT_ROOT / "01_bulk_vasp"
SLAB_VASP_DIR = PROJECT_ROOT / "04_slab_vasp"
METALS_VASP_DIR = PROJECT_ROOT / "05_metals_vasp"

# Aliases for backwards compatibility
SLAB_DIR = SLAB_QE_DIR
METALS_DIR = METALS_QE_DIR
VASP_BULK_VC_RELAX_DIR = BULK_VASP_DIR
VASP_SLAB_DIR = SLAB_VASP_DIR
VASP_METALS_DIR = METALS_VASP_DIR


def slab_profile_dir(profile_name):
    """Return the ``04_slab`` directory holding a named slab profile."""
    return SLAB_DIR / profile_name


def vasp_stage_root(stage="05_metals"):
    """Return the ``<stage>_vasp`` sibling directory (e.g. ``05_metals_vasp``)."""
    return PROJECT_ROOT / f"{stage}_vasp"


def vasp_metal_dir(metal, shape="Tetra"):
    """``05_metals_vasp/<Metal>/<Shape>`` for one cluster isomer."""
    return VASP_METALS_DIR / metal / shape


def vasp_adsorbed_dir(metal, shape, slab_profile, site):
    """``05_metals_vasp/<Metal>/<Shape>/Adsorbed/<slab>/<site>``."""
    return vasp_metal_dir(metal, shape) / "Adsorbed" / slab_profile / site
