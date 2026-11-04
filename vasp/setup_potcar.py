"""Assemble a VASP ``POTCAR`` from a PBE PAW library (pure stdlib, no ASE).

POTCARs are not stored in the repository (they are licensed VASP datasets,
exactly like the QE ``pseudo_dir`` lives only on the cluster).  Run this once
in each VASP work directory, pointed at your library, to materialise ``POTCAR``
from the per-species ``POTCAR_SPECIES`` file written by
:mod:`vasp.generate_inputs`.

Example::

    cd 05_metals_vasp/Co4/Tetra/Adsorbed/3_layers/lico
    python ../../../vasp/setup_potcar.py --pp-path /apps/pseudo/PBE
"""
import argparse
import sys
from pathlib import Path

DEFAULT_PP_PATH = None  # must be supplied or set via VASP_PP_PATH env


def assemble_potcar(species_path, out_path="POTCAR", pp_path=None, subdir=""):
    """Concatenate per-species POTCARs into ``out_path``.

    ``species_path`` points at the ``POTCAR_SPECIES`` file (one dataset name per
    line, e.g. ``Li Co O``).  Each dataset is read from ``pp_path/<subdir>/<name>/POTCAR``.
    """
    import os
    if pp_path is None:
        pp_path = os.environ.get("VASP_PP_PATH")
    if not pp_path:
        raise SystemExit("Set --pp-path or export VASP_PP_PATH to your PBE PAW library.")

    base = Path(pp_path)
    species = Path(species_path).read_text(encoding="utf-8", errors="ignore").split()
    chunks = []
    for name in species:
        candidate = base / subdir / name / "POTCAR" if subdir else base / name / "POTCAR"
        if not candidate.exists():
            # common vendor layout: <pp_path>/<FUNC>/<sym>/POTCAR
            candidate = base / "PBE" / name / "POTCAR"
        if not candidate.exists():
            raise FileNotFoundError(f"POTCAR for {name!r} not found under {base}")
        text = candidate.read_text(encoding="utf-8", errors="ignore").rstrip("\n") + "\n"
        chunks.append(text)
    Path(out_path).write_text("".join(chunks), encoding="utf-8")
    return out_path


def main(argv=None):
    ap = argparse.ArgumentParser(description="Assemble VASP POTCAR from a PBE PAW library.")
    ap.add_argument("--pp-path", default=None, help="root of your VASP PBE PAW library")
    ap.add_argument("--subdir", default="", help="optional sub-dir (e.g. 'PBE')")
    ap.add_argument("--species-file", default="POTCAR_SPECIES")
    ap.add_argument("--out", default="POTCAR")
    args = ap.parse_args(argv)
    result = assemble_potcar(args.species_file, args.out, args.pp_path, args.subdir)
    print(f"wrote {result}")


if __name__ == "__main__":
    main()
