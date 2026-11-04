"""Construction of VASP input decks.

This is the VASP counterpart to :mod:`lco.qe_input`: it builds the
``INCAR``, ``KPOINTS`` and ``POSCAR`` for a run from a structure plus a small
set of overrides, and guarantees the same reproducibility invariant that the
QE side enforces -- **every variant lives in its own directory with a unique
job name**, so VASP (``cwd``-based scratch) cannot silently clobber a sibling
variant the way QE ``prefix``/``outdir`` collisions used to.

Mapping summary (QE -> VASP), also documented in the README:

* ``vc-relax`` -> ``IBRION=1, ISIF=3`` (+ selective dynamics to freeze slab);
  ``relax`` -> ``IBRION=2, ISIF=2``; ``scf`` -> ``IBRION=-1, NSW=0, ISIF=2``.
* ``nspin=2`` -> ``ISPIN=2`` + ``LASPH=.TRUE.`` (+ ``MAGMOM`` in POSCAR).
* ``smearing='gaussian', degauss=0.01 Ry`` -> ``ISMEAR=0, SIGMA=0.136`` eV.
* ``vdw_corr='grimme-d3'`` -> ``IVDW=11`` (DFT-D3 zero-damping).
* ``K_POINTS 4 3 1 0 0 0`` -> VASP MP grid ``4 3 1`` shift ``0 0 0``.
* ``HUBBARD (ortho-atomic) U Co-3d 3.2`` -> ``LDAU=.TRUE., LDAUTYPE=2`` with
  ``LDAUL/LDAUU/LDAUJ`` emitted per species (see note in :mod:`lco.chemistry`).
"""
import re

from lco.chemistry import STARTING_MAGNETIZATION

# ---------------------------------------------------------------------------
# Element -> PBE PAW dataset name in your cluster library ($VASP_PP_PATH).
# These are the standard (no semi-core) PBE PAW datasets; swap an element to a
# ``_sv``/``_pv``/``_h`` variant by editing this table or passing `paw_dataset`.
# ---------------------------------------------------------------------------
POTCAR_DB = {
    "Li":  "Li",
    "Co":  "Co",
    "O":   "O",
    "Pt":  "Pt",
    "Ni":  "Ni",
    "Cu":  "Cu",
}

#: Default VASP INCAR values, overridden by the builder arguments per stage.
INCAR_DEFAULTS = {
    "SYSTEM": "LCO campaign (VASP PBE+PBE+U+D3)",
    "ISTART": 0,
    "ICHARG": 2,
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "EDIFFG": -1e-3,
    "ISMEAR": 0,
    "SIGMA": 0.136,
    "ISIF": 2,
    "IBRION": 2,
    "NSW": 0,
    "ISPIN": 2,
    "LASPH": True,
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": True,
    "IVDW": 11,
    "MLIFX": False,
    "MLFF": False,
    "MLU": False,
}

#: k-point divisions mirroring the QE campaign.
KPOINTS_SURFACE = (4, 3, 1)
KPOINTS_BULK = (12, 12, 4)

#: l-quantum integer for each angular-momentum letter (VASP LDAUL convention:
#: -1 none, 0 s, 1 p, 2 d, 3 f).
_L_QUANTUM = {"s": 0, "p": 1, "d": 2, "f": 3}


def _fmt(value):
    """Format a Python/INCAR value the way VASP expects."""
    if isinstance(value, bool):
        return ".TRUE." if value else ".FALSE."
    if isinstance(value, float):
        return repr(value)
    return str(value)


def variant_label(system, site, tag):
    """Unique, human-readable VASP job label for one experiment variant.

    Mirrors :func:`lco.qe_input.variant_prefix` (components: project marker,
    system tag, site, experiment tag).  VASP has no ``prefix`` keyword, so the
    unique directory + this ``job-name`` together enforce variant isolation --
    the VASP analogue of the QE prefix/outdir gotcha.
    """
    return f"lco_vasp_{system}_{site}_{tag}"


def potcar_species(elements):
    """Return the ordered list of PBE PAW dataset names for `elements`.

    The order matches the species order of the POSCAR so the concatenated
    ``POTCAR`` (built on the cluster from ``$VASP_PP_PATH``) lines up with
    ``POSCAR`` species, ``MAGMOM``, and ``LDAUL``/``LDAUU``/``LDAUJ``.
    """
    seen = []
    for el in elements:
        if el not in seen:
            seen.append(el)
    missing = [el for el in seen if el not in POTCAR_DB]
    if missing:
        raise KeyError(f"Element(s) {missing} missing in POTCAR_DB configuration.")
    return [POTCAR_DB[el] for el in seen]


def _ldaU_incar(elements, ldau=None):
    """Return the LDAU/LDAUL/LDAUU/LDAUJ lines for the given species order.

    `ldau` maps element -> U (eV).  Elements absent from `ldau` get LDAUL=-1
    (no on-site term).  LDAUJ defaults to 0.0 (U_eff = U - J = LDAUU).
    """
    if not ldau:
        return {}
    order = []
    for el in elements:
        if el not in order:
            order.append(el)
    ldaul, ldauu, ldauj = [], [], []
    for el in order:
        u = ldau.get(el, 0.0)
        ldaul.append(str(_L_QUANTUM["d"] if u else -1))   # U applied to d only
        ldauu.append(_fmt(u))
        ldauj.append(_fmt(0.0))
    return {
        "LDAU": True,
        "LDAUTYPE": 2,
        "LDAUL": " ".join(ldaul),
        "LDAUU": " ".join(ldauu),
        "LDAUJ": " ".join(ldauj),
    }


def build_incar(elements, *, calculation="relax", encut=520, sigma=0.136,
                ldau=None, overrides=None, **kw):
    """Build a VASP ``INCAR`` string.

    Parameters
    ----------
    elements
        ordered species list (defines POTCAR/LDAU order).
    calculation
        ``'vc-relax'`` -> ISIF=3, IBRION=1, NSW>=1;
        ``'relax'``    -> ISIF=2, IBRION=2;
        ``'scf'``      -> IBRION=-1, NSW=0.
    ldau
        dict element -> U(eV).  Pass ``HUBBARD_U`` (from :mod:`lco.chemistry`)
        to reproduce the QE campaign's Co U=3.2.
    overrides
        extra INCAR keys merged after the defaults (last word wins).
    **kw
        individual INCAR overrides (e.g. ``ISIF=3``).
    """
    cfg = dict(INCAR_DEFAULTS)
    if calculation == "vc-relax":
        cfg.update(IBRION=1, ISIF=3, NSW=kw.pop("NSW", 100))
    elif calculation == "relax":
        cfg.update(IBRION=2, ISIF=2, NSW=kw.pop("NSW", 200))
    elif calculation == "scf":
        cfg.update(IBRION=-1, ISIF=2, NSW=0)
    else:
        raise ValueError(f"Unknown VASP calculation='{calculation}'.")

    cfg.update(ENCUT=encut, SIGMA=sigma)
    cfg.update(kw)
    if ldau:
        cfg.update(_ldaU_incar(elements, ldau))
    if overrides:
        cfg.update(overrides)

    lines = [f"SYSTEM = {cfg.pop('SYSTEM', INCAR_DEFAULTS['SYSTEM'])}"]
    for key, value in cfg.items():
        lines.append(f"{key} = {_fmt(value)}")
    lines.append("")
    return "\n".join(lines) + "\n"


def build_kpoints(grid, shift=(0, 0, 0), comment=None):
    """Build a Monkhorst-Pack ``KPOINTS`` file (automatic mesh)."""
    g = " ".join(str(int(n)) for n in grid)
    s = " ".join(str(int(n)) for n in shift)
    header = comment or f"Automatic mesh {g}, shift {s}"
    return f"{header}\n0\n{g}\n{s}\n"


def build_kpoints_gamma(comment="Gamma-only"):
    """Build a Gamma-only ``KPOINTS`` file (isolated cluster / molecule)."""
    return f"{comment.strip()}\n0\nGamma\n1 1 1\n0 0 0\n"


def build_poscar(atoms, comment=None, direct=True, selective_dynamics=None,
                 masses=None):
    """Build a VASP ``POSCAR`` string from an ASE ``Atoms`` (or duck-typed).

    `atoms` must provide ``get_cell()``, ``get_scaled_positions()``,
    ``get_chemical_symbols()`` and ``get_masses()`` (the ASE Atoms API), which
    keeps this function testable without ASE and independent of ASE's own
    (version-dependent) VASP writer.
    """
    cell = atoms.get_cell()
    symbols = list(atoms.get_chemical_symbols())

    # unique element order for the POSCAR species line
    species, counts = [], []
    for s in symbols:
        if s in species:
            counts[species.index(s)] += 1
        else:
            species.append(s)
            counts.append(1)

    lines = []
    lines.append(comment or " ".join(symbols))
    lines.append("   1.0")
    for vec in cell:
        lines.append(f"  {vec[0]:21.16f}  {vec[1]:21.16f}  {vec[2]:21.16f}")
    lines.append("  " + "  ".join(species))
    lines.append("  " + "  ".join(str(c) for c in counts))
    coords = atoms.get_scaled_positions() if direct else atoms.get_positions()
    use_seldyn = selective_dynamics is not None and len(selective_dynamics) == len(symbols)
    if use_seldyn:
        lines.append("Selective dynamics")
    lines.append("direct" if direct else "cartesian")
    if use_seldyn:
        flag = lambda fix: ("F" if fix else "T")
        for (x, y, z), fix in zip(coords, selective_dynamics):
            lines.append(f"  {x:21.16f} {y:21.16f} {z:21.16f}  {flag(fix)} {flag(fix)} {flag(fix)}")
    else:
        for x, y, z in coords:
            lines.append(f"  {x:21.16f} {y:21.16f} {z:21.16f}")
    lines.append("")
    return "\n".join(lines)


def format_magmom(atoms, ldau_elements=None):
    """Return ``(magmom_str, magmom_list)`` sized to the atom count for INCAR.

    Default guess: the value in :data:`lco.chemistry.STARTING_MAGNETIZATION`
    for each element (Co ~0.6 uB), 0.0 elsewhere.  Seeds override per-system
    to reproduce the QE ``starting_magnetization`` values.
    """
    symbols = list(atoms.get_chemical_symbols())
    # MAGMOM is tied to ISPIN, not to LDAU: every atom gets a guess from the
    # table (0.0 for non-magnetic elements).  `ldau_elements` is accepted for
    # API symmetry but does not change the result.
    mag = [STARTING_MAGNETIZATION.get(s, 0.0) for s in symbols]
    return " ".join(f"{m:.2f}" for m in mag), mag
