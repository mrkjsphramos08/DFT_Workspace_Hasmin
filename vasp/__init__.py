"""VASP input-generation / job-submission helpers.

Mirrors the top-level Quantum ESPRESSO generators (``Adsorbed_QE_*.py``,
``run_sh*.py``) but for VASP:

* :mod:`vasp.generate_inputs` builds ``INCAR``/``KPOINTS``/``POSCAR`` per site.
* :mod:`vasp.run_sh` emits SLURM scripts that run ``vasp_std``.

Both reuse the shared :mod:`lco` library so QE and VASP stay in lock-step.
"""
