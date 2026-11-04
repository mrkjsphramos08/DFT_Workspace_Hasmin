"""Engine-agnostic chemical reference data for the LCO campaign.

This module is the single source of truth for the chemical constants that are
shared between the Quantum ESPRESSO (QE) and VASP stages, so the two engines
never drift apart:

* :data:`ELEMENT_MASSES`  -- atomic masses used when writing ``ATOMIC_SPECIES``
  (QE) or the ``POTCAR``-free mass line / POSCAR validation (VASP).
* :data:`HUBBARD_U`       -- the campaign's on-site Coulomb ``U`` (eV) per
  element.  In QE this is applied as ``HUBBARD (ortho-atomic)`` (see
  :data:`lco.qe_input.HUBBARD_BLOCK`); in VASP it is mapped onto the
  Dudarev parameters ``LDAUU``/``LDAUJ`` (:mod:`lco.vasp_input`).

Note on the U convention
------------------------
The QE value ``U Co-3d 3.2`` is an *ortho-atomic* (atomic-limit) U.  VASP's
Dudarev ``LDAUU`` is ``U_eff = U - J``.  Setting ``LDAUJ = 0`` makes
``U_eff = LDAUU = 3.2`` eV, which is the closest mechanical analogue; the two
are *not* guaranteed to be numerically identical and the VASP runs should be
validated against the QE reference.  ``LDAUJ`` can be raised to a finite
 screening value (e.g. 0.7-0.9 eV for Co) if you prefer a non-zero J.
"""
# Atomic masses (amu) -- kept in sync with lco.qe_input.PSEUDO_DB.
ELEMENT_MASSES = {
    "Li": 6.941,
    "Co": 58.933,
    "O": 15.999,
    "Pt": 195.084,
    "Ni": 58.693,
    "Cu": 63.546,
}

#: Hubbard U (eV) applied per element.  Only Co is treated here, matching the
#: ``HUBBARD (ortho-atomic) U Co-3d 3.2`` card used throughout the QE inputs.
HUBBARD_U = {
    "Co": 3.2,
}

#: Starting magnetic moments (uB) per element for VASP ``MAGMOM``.  These are
#: deliberately *guesses* matching the magnitude of the QE campaign's
#: ``starting_magnetization`` tags (Co ~0.1 in bulk LCO, ~0.8 on the bare
#: Co4 cluster).  Seeds below override these per-system to reproduce the QE
#: values exactly; the generator uses this table as the default.
STARTING_MAGNETIZATION = {
    "Li": 0.0, "Co": 0.6, "O": 0.0, "Pt": 1.0, "Ni": 0.6, "Cu": 0.6,
}
