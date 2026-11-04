"""Shared, ASE-based structure I/O for the QE and VASP stages.

ASE is imported **lazily** inside every function (the import is cached by the
interpreter after the first call), so the module can be imported, byte-compiled
and unit-tested for its pure-Python helpers in environments that do not have
ASE installed -- mirroring the pattern already used by :mod:`lco.qe_output`.
"""
from pathlib import Path


def _ase():
    """Return ASE's read/write pair, importing ASE on first use."""
    from ase.io import read as ase_read, write as ase_write
    return ase_read, ase_write


def read_structure(path):
    """Read any structure format ASE understands (cif, vasp, xyz, ...) into Atoms."""
    ase_read, _ = _ase()
    return ase_read(str(path))


def write_structure(atoms, path, fmt=None):
    """Write an ``Atoms`` object, guessing the format from the extension."""
    _, ase_write = _ase()
    name = Path(path).suffix.lower()
    if fmt is None:
        fmt = {"cif": "cif", "vasp": "vasp", "poscar": "vasp", "xyz": "xyz"}.get(name.lstrip("."), "cif")
    ase_write(str(path), atoms, format=fmt)


def freeze_by_z(atoms, z_cutoff_fraction=0.5, num_above=0):
    """Return a per-atom selective-dynamics mask.

    ``True`` means the atom is *fixed* (frozen).  Mirrors the QE workflow's
    "freeze the bottom half of the slab" rule: atoms in the bottom
    ``z_cutoff_fraction`` of the mapped slab height are fixed, except the
    top ``num_above`` atom layers (cluster atoms sitting above the slab) which
    are always left free.
    """
    import numpy as np
    n = len(atoms)
    n_slab = n - num_above
    scaled_z = np.asarray(atoms.get_scaled_positions())[:, 2]
    mapped_z = (scaled_z + 0.1) % 1.0 - 0.1  # handle wrap-around near 0/1
    slab_z = mapped_z[:n_slab]
    z_min, z_max = slab_z.min(), slab_z.max()
    z_thr = z_min + z_cutoff_fraction * (z_max - z_min)

    fixed = [False] * n
    for i in range(n_slab):
        if mapped_z[i] < z_thr:
            fixed[i] = True
    return fixed
