"""Parsing of VASP output files.

Following the same "parse the text directly, do not trust a library reader"
rule the project already uses for Quantum ESPRESSO (:mod:`lco.qe_output`), this
module extracts the few quantities the campaign needs from the VASP text
outputs (``OSZIC`` / ``OUTCAR``):

* the total/free energy per electronic step (``TOTEN``);
* the cell magnetization;
* whether electronic convergence was reached.

Everything here is plain string parsing with the standard library, so the
module imports (and its unit tests) need neither ASE nor VASP.
"""
import re
from pathlib import Path

# --------------------------------------------------------------------------
# Regexes -- written against the real VASP OUTCAR/OSZIC formatting, but tolerant
# of the extra whitespace VASP is fond of.
# --------------------------------------------------------------------------

#   free  energy   TOTEN                 =    -25.1234567      eV
_ENERGY_RE = re.compile(r"free\s+energy\s+TOTEN\s*=\s*(-?\d+\.\d+)")

#   total magnetization =    3.000  Bohr mag/cell
_MAGNET_RE = re.compile(r"total magnetization\s*=\s*(-?[\d.]+)\s*Bohr", re.IGNORECASE)

# VASP writes this line when an SCF cycle hits the requested accuracy.
_CONVERGED_RE = re.compile(r"reached required accuracy", re.IGNORECASE)

#   reached required accuracy - stopping electronic iterations
_ELEC_CONVERGED_RE = re.compile(
    r"reached required accuracy - stopping electronic iterations", re.IGNORECASE)

# An SCF step that did not converge within electron max steps:
#   "TOO MANY BANS" / "reached required accuracy" absent but "General timing" present
_JOB_DONE_RE = re.compile(r"General timing and performance information")


def _read_text(path):
    return Path(path).read_text(encoding="utf-8", errors="ignore")


def parse_total_energies(path):
    """Return all ``free energy TOTEN`` values (one per electronic step)."""
    text = _read_text(path)
    return [float(m) for m in _ENERGY_RE.findall(text)]


def parse_final_energy(path):
    """Return the last recorded free energy (eV), or ``None``."""
    energies = parse_total_energies(path)
    return energies[-1] if energies else None


def parse_total_magnetization(path):
    """Return the last cell magnetization (Bohr mag/cell), or ``None``."""
    text = _read_text(path)
    mags = _MAGNET_RE.findall(text)
    return float(mags[-1]) if mags else None


def parse_convergence(path):
    """Return the convergence summary of a VASP run as a dict.

    Keys (mirroring :func:`lco.qe_output.parse_scf_history`):

    ``energies``
        free-energy TOTEN at every electronic step (eV);
    ``final_energy``
        the last energy (eV), or ``None``;
    ``magnetization``
        last cell magnetization (Bohr mag/cell), or ``None``;
    ``converged``
        whether electronic convergence was reached;
    ``iterations``
        number of electronic steps recorded;
    ``job_done``
        whether the run reached the end-of-run timing banner.
    """
    text = _read_text(path)
    energies = [float(m) for m in _ENERGY_RE.findall(text)]
    mags = _MAGNET_RE.findall(text)
    return {
        "path": str(path),
        "energies": energies,
        "iterations": len(energies),
        "final_energy": energies[-1] if energies else None,
        "magnetization": [float(m) for m in mags],
        "converged": bool(_CONVERGED_RE.search(text)),
        "job_done": bool(_JOB_DONE_RE.search(text)),
    }


def final_energy(history):
    """Return the best (last) energy from a :func:`parse_convergence` dict."""
    energies = history.get("energies") or []
    return energies[-1] if energies else None
