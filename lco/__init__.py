"""Shared library for the LCO bulk adsorption / Quantum ESPRESSO campaign.

The package collects the pieces that were previously duplicated across the
top-level scripts:

``lco.paths``
    Single source of truth for repository path resolution.
``lco.qe_output``
    Parsing of Quantum ESPRESSO outputs (final coordinates, SCF history).
``lco.qe_input``
    Building of ``pw.x`` input decks from a template plus parameter overrides.

Modules ``lco.sites``, ``lco.slurm``, ``lco.variants`` and ``lco.convergence``
are planned but not yet implemented.
"""

__all__ = [
    "paths",
    "qe_output",
    "qe_input",
]
