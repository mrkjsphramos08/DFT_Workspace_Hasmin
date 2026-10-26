"""Shared library for the LCO bulk adsorption / Quantum ESPRESSO campaign.

The package collects the pieces that were previously duplicated across the
top-level scripts:

``lco.paths``
    Single source of truth for repository path resolution.
``lco.sites``
    Canonical adsorption-site definitions and the generic site resolver.
``lco.qe_output``
    Parsing of Quantum ESPRESSO outputs (final coordinates, SCF history).
``lco.qe_input``
    Building of ``pw.x`` input decks from a template plus parameter overrides.
``lco.slurm``
    Generation of SLURM submission scripts.
``lco.variants``
    Named convergence-experiment parameter variants.
``lco.convergence``
    Diagnosis of SCF convergence failures and next-step recommendations.
"""

__all__ = [
    "paths",
    "sites",
    "qe_output",
    "qe_input",
    "slurm",
    "variants",
    "convergence",
]
