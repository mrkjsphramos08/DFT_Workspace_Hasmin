"""Parsing of Quantum ESPRESSO output files.

Two things are needed from a ``pw.x`` output:

* the **final relaxed structure**, used to rebuild CIFs.  ASE's
  ``espresso-out`` reader fails on a number of the outputs in this repository
  (it raises ``AssertionError`` on the ``3_layers_new`` slab), so the
  ``Begin final coordinates`` block is parsed directly.  This is the reliable
  route and removes the need for the ``espresso-xml`` fallback that never
  worked -- ASE has no ``espresso-xml`` reader at all.
* the **SCF history**, used to diagnose why a run failed to converge.
"""
import re

# --------------------------------------------------------------------------
# Final coordinates
# --------------------------------------------------------------------------

_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?"
_UNIT = r"[\w.]+"
_FINAL_COORDS_RE = re.compile(
    rf"^\s*ATOMIC_POSITIONS\s*(?:\(\s*({_UNIT})\s*\)|({_UNIT}))", re.IGNORECASE
)
_CELL_RE = re.compile(
    rf"^\s*CELL_PARAMETERS\s*(?:\(\s*({_UNIT})\s*\)|({_UNIT}))",
    re.IGNORECASE | re.MULTILINE,
)
_ELEM_RE = re.compile(rf"^([A-Z][a-z]?)\s+({_NUM})\s+({_NUM})\s+({_NUM})")


def _fortran_float(text):
    """Parse a Fortran-style float (``1.0D-4`` or ``1.0E-4``)."""
    return float(text.replace("D", "E").replace("d", "e"))


def _fractional_to_cartesian(positions, cell):
    """Convert scaled (``crystal``) positions to Cartesian using the cell rows."""
    import numpy as np
    frac = np.asarray(positions, dtype=float)
    lattice = np.asarray(cell, dtype=float)
    return [tuple(row) for row in frac @ lattice]


def _to_angstrom(values):
    from ase.units import Bohr
    return [v * Bohr for v in values]


def _identity(values):
    return list(values)


def _parse_vector_block(lines, converter):
    """Parse a 3-row numeric block into a list of 3-tuples.

    Leading blank lines are skipped (the block normally starts on the line
    after its header) and parsing stops at the first line that is not three
    numbers.
    """
    rows = []
    for line in lines:
        if not rows and not line.strip():
            continue  # skip blank lines before the block starts
        parts = line.split()
        if len(parts) < 3:
            break
        try:
            rows.append(tuple(_fortran_float(p) for p in parts[:3]))
        except ValueError:
            break
    if len(rows) < 3:
        return None
    return [converter(list(v)) for v in rows[:3]]

def read_final_coordinates(path, cell=None):
    """Return the final structure of a ``pw.x`` output as an ASE ``Atoms``.

    Parameters
    ----------
    path
        the ``pw.x`` output file.
    cell
        optional cell (list of three 3-tuples, Angstrom).  When given it is
        used verbatim; this is the preferred route because the cell in the
        output header is only printed to limited precision.

    Returns ``None`` when the output contains no ``Begin final coordinates``
    block (for example when the job died before finishing the first SCF).
    """
    from ase import Atoms

    text = _read_text(path)
    if "Begin final coordinates" not in text:
        return None

    tail = text.split("Begin final coordinates", 1)[1]
    lines = tail.splitlines()

    # Find the ATOMIC_POSITIONS / CELL_PARAMETERS headers and their units.
    positions_unit = "angstrom"
    cell_unit = "angstrom"
    pos_idx = cell_idx = None
    for i, line in enumerate(lines):
        m = _FINAL_COORDS_RE.match(line)
        if m and pos_idx is None:
            positions_unit = m.group(1) or m.group(2)
            pos_idx = i
        m = _CELL_RE.match(line)
        if m and cell_idx is None:
            cell_unit = m.group(1) or m.group(2)
            cell_idx = i

    if pos_idx is None:
        return None

    pos_conv = _to_angstrom if positions_unit in ("a.u.", "bohr") else _identity
    cell_conv = _to_angstrom if cell_unit in ("a.u.", "bohr") else _identity
    positions_fractional = positions_unit.lower() == "crystal"

    # ATOMIC_POSITIONS: read until a line that is not an element + 3 numbers.
    symbols, positions = [], []
    for line in lines[pos_idx + 1:]:
        if not line.strip():
            if symbols:
                break
            continue
        m = _ELEM_RE.match(line.strip())
        if not m:
            break
        symbols.append(m.group(1))
        positions.append(tuple(_fortran_float(m.group(i)) for i in (2, 3, 4)))

    if not symbols:
        return None

    if cell is None:
        # The cell may be reported inside the final-coordinates block ...
        if cell_idx is not None:
            cell = _parse_vector_block(lines[cell_idx + 1:], cell_conv)
        if cell is None:
            # ... otherwise fall back to the header information.
            cell = _cell_from_crystal_axes(text, cell_conv)
        if cell is None:
            cell = _cell_from_input_echo(text, cell_conv)

    if cell is None:
        return None

    if positions_fractional:
        positions = _fractional_to_cartesian(positions, cell)

    atoms = Atoms(symbols=symbols, positions=pos_conv(positions), cell=cell, pbc=True)
    return atoms


def _cell_from_input_echo(text, converter):
    """Recover the cell from the ``CELL_PARAMETERS`` block QE echoes in the output."""
    matches = list(re.finditer(r"CELL_PARAMETERS\s*\(\s*(\w+)\s*\)", text, re.IGNORECASE))
    if not matches:
        return None
    last = matches[-1]
    lines = text[last.end():].splitlines()
    block = []
    for line in lines:
        parts = line.split()
        if len(parts) < 3:
            break
        try:
            block.append(tuple(_fortran_float(p) for p in parts[:3]))
        except ValueError:
            break
    if len(block) < 3:
        return None
    return [converter(list(v)) for v in block[:3]]

def _cell_from_crystal_axes(text, converter):
    """Rebuild the cell from QE's ``crystal axes`` / ``lattice parameter`` header.

    A fixed-cell ``relax`` does not print a ``CELL_PARAMETERS`` block, so the
    only cell information in the output is the crystal axes given in units of
    ``alat``.  The printed precision limits the result to roughly 1e-5 A, which
    is why :func:`read_cell_from_input` is preferred when the input deck is
    available.
    """
    import numpy as np
    from ase.units import Bohr

    alat_match = re.search(r"lattice parameter \(alat\)\s*=\s*([\d.]+)", text)
    if not alat_match:
        return None
    alat = float(alat_match.group(1)) * Bohr  # Angstrom

    axes = []
    for match in re.finditer(r"a\((\d)\)\s*=\s*\(\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*\)", text):
        axes.append([float(match.group(i)) for i in (2, 3, 4)])
        if len(axes) == 3:
            break
    if len(axes) != 3:
        return None

    cell = np.array(axes) * alat
    return [converter(list(row)) for row in cell]


def read_cell_from_input(input_path):
    """Read the ``CELL_PARAMETERS`` block of a ``pw.x`` input deck.

    Returns the cell in Angstrom as a list of three 3-tuples, or ``None`` when
    the deck does not declare one (for example when it uses ``ibrav``).
    """
    text = _read_text(input_path)
    match = _CELL_RE.search(text)
    if not match:
        return None
    unit = match.group(1) or match.group(2)
    converter = _to_angstrom if unit in ("a.u.", "bohr") else _identity
    return _parse_vector_block(text[match.end():].splitlines(), converter)

# --------------------------------------------------------------------------
# SCF history
# --------------------------------------------------------------------------

_ACCURACY_RE = re.compile(r"estimated scf accuracy\s*<\s*([\d.eE+-]+)")
_CONVERGED_RE = re.compile(r"convergence has been achieved in\s+(\d+)\s+iterations")
_NOT_CONVERGED_RE = re.compile(r"convergence NOT achieved after\s+(\d+)\s+iterations")
_THRESHOLD_RE = re.compile(r"scf convergence threshold\s+=\s*(\S+)")
_BETA_RE = re.compile(r"mixing beta\s*=\s*([\d.]+)")
_NDIM_RE = re.compile(r"number of iterations used\s*=\s*(\d+)\s+(\S+)\s+mixing")
_MAXSTEP_RE = re.compile(r"max number of scf iterations\s*=\s*(\d+)")
_MIXING_MODE_RE = re.compile(r"Mixing scheme\s*=\s*(\S+)")
_NEGATIVE_RHO_RE = re.compile(r"negative rho")
_MAGN_RE = re.compile(r"total magnetization\s+=\s+(-?[\d.]+)\s+Bohr mag/cell")
_JOB_DONE_RE = re.compile(r"JOB DONE")
_ENERGY_RE = re.compile(r"^\s*!\s*total energy\s+=\s+(-?[\d]+\.\d+)\s+Ry", re.MULTILINE)


def _read_text(path):
    from pathlib import Path
    return Path(path).read_text(encoding="utf-8", errors="ignore")


def parse_scf_history(path):
    """Return the SCF convergence history of a ``pw.x`` output.

    The returned dict carries:

    ``accuracies``
        the ``estimated scf accuracy`` after every SCF iteration;
    ``converged``
        whether QE reported successful convergence;
    ``threshold``
        the ``scf convergence threshold`` that was in force;
    ``beta`` / ``mixing_ndim`` / ``mixing_mode``
        the mixing parameters actually used;
    ``max_steps``
        ``electron_maxstep``;
    ``iterations``
        number of SCF iterations performed;
    ``negative_rho``
        whether the output shows a negative charge density;
    ``job_done``
        whether the run reached the ``JOB DONE`` banner.
    """
    text = _read_text(path)

    accuracies = [float(m) for m in _ACCURACY_RE.findall(text)]
    converged = _CONVERGED_RE.search(text)
    not_converged = _NOT_CONVERGED_RE.search(text)
    threshold = _THRESHOLD_RE.search(text)
    beta = _BETA_RE.search(text)
    ndim = _NDIM_RE.search(text)
    max_steps = _MAXSTEP_RE.search(text)
    mixing_mode = _MIXING_MODE_RE.search(text)

    energy = _ENERGY_RE.findall(text)

    return {
        "path": str(path),
        "accuracies": accuracies,
        "iterations": len(accuracies),
        "converged": bool(converged) and not not_converged,
        "iterations_to_converge": int(converged.group(1)) if converged else None,
        "threshold": float(threshold.group(1)) if threshold else None,
        "beta": float(beta.group(1)) if beta else None,
        "mixing_ndim": int(ndim.group(1)) if ndim else None,
        "mixing_mode": mixing_mode.group(1) if mixing_mode else None,
        "max_steps": int(max_steps.group(1)) if max_steps else None,
        "negative_rho": bool(_NEGATIVE_RHO_RE.search(text)),
        "job_done": bool(_JOB_DONE_RE.search(text)),
        "final_energy": float(energy[-1]) if energy else None,
        "magnetization": [float(m) for m in _MAGN_RE.findall(text)],
    }


def final_accuracy(history):
    """Return the best (lowest) SCF accuracy reached, or ``None``."""
    acc = history.get("accuracies") or []
    return min(acc) if acc else None
