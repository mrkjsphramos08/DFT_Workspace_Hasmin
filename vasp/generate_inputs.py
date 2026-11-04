"""Generate VASP input decks (``INCAR``/``KPOINTS``/``POSCAR``) from CIFs.

This is the VASP analogue of ``Adsorbed_QE_lico.py`` / ``Adsorbed_QE_Pt.py``:
it reads a *relaxed* structure (a CIF produced by the QE side) and writes a
VASP deck into the parallel ``<stage>_vasp/`` tree.  Reusing the relaxed QE
structures keeps the two engines locked to the same geometry; only the
electronic-structure method changes.

Usage (adsorbed campaign)::

    python -m vasp.generate_inputs --metal Co4 --shape Tetra \\
        --slab 3_layers --site lico --calculation relax --ldaU Co=3.2

    # everything:
    python -m vasp.generate_inputs --all

    # one arbitrary structure:
    python -m vasp.generate_inputs --cif path/to/struc.cif --out-dir work/vasp_struc
"""
import argparse
import argparse
from pathlib import Path

from lco.paths import METALS_DIR, VASP_METALS_DIR
from lco.ase_io import read_structure
from lco.vasp_input import (
    build_incar, build_kpoints, build_poscar, potcar_species,
    format_magmom, KPOINTS_SURFACE, KPOINTS_BULK,
)
from adsorption_sites import VALID_SITES
from adsorption_config import SLAB_CHOICES


def _parse_ldau(spec):
    """Parse ``Co=3.2,Ni=6.0`` -> {'Co': 3.2, 'Ni': 6.0}.  Empty/None -> None."""
    if not spec:
        return None
    out = {}
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        el, val = token.split("=", 1)
        out[el.strip()] = float(val)
    return out


def write_vasp_inputs(out_dir, *, atoms=None, cif_path=None, calculation="relax",
                      kpoints=KPOINTS_SURFACE, encut=520, sigma=0.136, ldau=None,
                      selective_dynamics=None, comment=None, magmom=None,
                      overwrite=False):
    """Write INCAR/KPOINTS/POSCAR (+ a .cif copy) into ``out_dir``.

    Parameters
    ----------
    atoms
        an ASE Atoms object.  Pass this to skip ASE and unit-test without ASE.
    cif_path
        source CIF; read with ASE only when ``atoms`` is None.
    selective_dynamics
        per-atom booleans (True = fixed), e.g. from ``lco.ase_io.freeze_by_z``.
    magmom
        per-atom starting magnetizations; emitted as VASP ``MAGMOM``.
    """
    out_dir = Path(out_dir)
    done = out_dir / "INCAR"
    if done.exists() and not overwrite:
        return [done]
    out_dir.mkdir(parents=True, exist_ok=True)

    if atoms is None:
        if cif_path is None:
            raise ValueError("Either `atoms` or `cif_path` must be provided.")
        atoms = read_structure(cif_path)

    symbols = list(atoms.get_chemical_symbols())

    extra = {}
    if magmom:
        extra["MAGMOM"] = " ".join(f"{m:.2f}" for m in magmom)
    incar = build_incar(symbols, calculation=calculation, encut=encut,
                        sigma=sigma, ldau=ldau, **extra)
    kpoints_str = build_kpoints(kpoints)
    poscar = build_poscar(atoms, comment=comment or " ".join(symbols),
                          selective_dynamics=selective_dynamics)

    (out_dir / "INCAR").write_text(incar, encoding="utf-8")
    (out_dir / "KPOINTS").write_text(kpoints_str, encoding="utf-8")
    (out_dir / "POSCAR").write_text(poscar, encoding="utf-8")
    (out_dir / "POTCAR_SPECIES").write_text(" ".join(potcar_species(symbols)) + "\n",
                                            encoding="utf-8")

    # provenance: keep the source structure alongside the VASP inputs.
    if cif_path:
        src = Path(cif_path).name
        (out_dir / src).write_text(Path(cif_path).read_text(encoding="utf-8",
                                                           errors="ignore"),
                                   encoding="utf-8")
    return [out_dir / "INCAR", out_dir / "KPOINTS", out_dir / "POSCAR",
            out_dir / "POTCAR_SPECIES"]


def _qe_metal_source(metals_dir, metal, shape, slab, site):
    """Find the QE source CIF for one adsorbed site."""
    site_dir = metals_dir / metal / shape / "Adsorbed" / slab / site
    cifs = sorted(site_dir.glob("adsorbed_*.cif"))
    if not cifs:
        raise FileNotFoundError(f"No adsorbed_*.cif in {site_dir}")
    return cifs[0]


def generate_adsorbed(metals_dir, vasp_dir, slab, sites, calculation, ldau,
                      kpoints, encut, sigma, freeze_z, overwrite, verbose=True):
    """Fan out VASP inputs over <Metal>/<Shape>/Adsorbed/<slab>/<site>/."""
    from lco.ase_io import freeze_by_z
    count = 0
    for metal_dir in sorted(metals_dir.iterdir()):
        if not metal_dir.is_dir():
            continue
        for shape_dir in sorted(m for m in metal_dir.iterdir() if m.is_dir()):
            ads = shape_dir / "Adsorbed" / slab
            if not ads.is_dir():
                continue
            for site_dir in sorted(s for s in ads.iterdir() if s.is_dir()):
                if site_dir.name not in sites:
                    continue
                try:
                    src = _qe_metal_source(metals_dir, metal_dir.name, shape_dir.name, slab, site_dir.name)
                except FileNotFoundError:
                    continue
                rel = Path(metal_dir.name) / shape_dir.name / "Adsorbed" / slab / site_dir.name
                out = vasp_dir / rel
                atoms = read_structure(src)
                sd = freeze_by_z(atoms, z_cutoff_fraction=freeze_z) if freeze_z else None
                _, mag = format_magmom(atoms)
                write_vasp_inputs(out, atoms=atoms, calculation=calculation,
                                  kpoints=kpoints, encut=encut, sigma=sigma,
                                  ldau=ldau, selective_dynamics=sd,
                                  comment=f"{metal_dir.name}_{shape_dir.name}@{site_dir.name}",
                                  magmom=mag, overwrite=overwrite)
                count += 1
                if verbose:
                    print(f"  wrote {out}")
    return count


def generate_clusters(metals_dir, vasp_dir, calculation, ldau, kpoints,
                      encut, sigma, overwrite, verbose=True):
    """Fan out VASP inputs over the bare <Metal>/<Shape> clusters."""
    count = 0
    for metal_dir in sorted(metals_dir.iterdir()):
        if not metal_dir.is_dir():
            continue
        for shape_dir in sorted(m for m in metal_dir.iterdir() if m.is_dir()):
            cifs = sorted(shape_dir.glob("*_latest.cif"))
            if not cifs:
                continue
            out = vasp_dir / metal_dir.name / shape_dir.name
            atoms = read_structure(cifs[0])
            _, mag = format_magmom(atoms)
            write_vasp_inputs(out, atoms=atoms, cif_path=cifs[0], calculation=calculation,
                              kpoints=kpoints, encut=encut, sigma=sigma,
                              ldau=ldau, comment=f"{metal_dir.name}_{shape_dir.name}",
                              magmom=mag, overwrite=overwrite)
            count += 1
            if verbose:
                print(f"  wrote {out}")
    return count


def main(argv=None):
    ap = argparse.ArgumentParser(description="Generate VASP input decks for the LCO campaign.")
    ap.add_argument("--cif", help="single source CIF to convert")
    ap.add_argument("--out-dir", help="output directory (with --cif)")
    ap.add_argument("--all", action="store_true", help="generate everything (clusters + adsorbed)")
    ap.add_argument("--metals-dir", default=str(METALS_DIR))
    ap.add_argument("--vasp-dir", default=str(VASP_METALS_DIR))
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--adsorbed", action="store_true", help="only adsorbed sites")
    mode.add_argument("--clusters", action="store_true", help="only bare clusters")
    ap.add_argument("--metal", help="filter: e.g. Co4 / Pt4 / Pt1")
    ap.add_argument("--shape", choices=["Tetra", "Rhom", "Rhom2"], default=None)
    ap.add_argument("--slab", choices=SLAB_CHOICES, default="3_layers")
    ap.add_argument("--site", choices=list(VALID_SITES), default=None)
    ap.add_argument("--calculation", choices=["vc-relax", "relax", "scf"], default="relax")
    ap.add_argument("--kpoints", type=int, nargs="+", default=None, help="e.g. 4 3 1")
    ap.add_argument("--ldaU", default="Co=3.2", help="e.g. Co=3.2  (use 'off' for none)")
    ap.add_argument("--encut", type=float, default=520)
    ap.add_argument("--sigma", type=float, default=0.136)
    ap.add_argument("--freeze-below", type=float, default=0.5, help="z-fraction to freeze (slab bottom)")
    ap.add_argument("--no-freeze", action="store_true", help="do not write selective dynamics")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args(argv)

    ldau = None if args.ldau == "off" else _parse_ldau(args.ldau)
    kpoints = (tuple(args.kpoints) if args.kpoints
               else KPOINTS_BULK if args.calculation == "vc-relax"
               else KPOINTS_SURFACE)
    metals_dir = Path(args.metals_dir)
    vasp_dir = Path(args.vasp_dir)

    if args.cif:
        out = Path(args.out_dir) if args.out_dir else vasp_dir / "single"
        written = write_vasp_inputs(out, cif_path=args.cif, calculation=args.calculation,
                                    kpoints=kpoints, encut=args.encut, sigma=args.sigma,
                                    ldau=ldau, overwrite=args.overwrite)
        print(f"wrote {len(written)} files -> {out}")
        return

    freeze_z = 0.0 if args.no_freeze else args.freeze_below
    if args.all or (not args.adsorbed and not args.clusters and not args.metal and not args.site):
        n = generate_clusters(metals_dir, vasp_dir, args.calculation, ldau, kpoints,
                              args.encut, args.sigma, args.overwrite)
        n += generate_adsorbed(metals_dir, vasp_dir, args.slab, set(VALID_SITES),
                               args.calculation, ldau, kpoints, args.encut,
                               args.sigma, freeze_z, args.overwrite)
        print(f"generated {n} VASP input sets")
        return

    if args.clusters:
        n = generate_clusters(metals_dir, vasp_dir, args.calculation, ldau, kpoints,
                              args.encut, args.sigma, args.overwrite)
    else:
        sites = {args.site} if args.site else set(VALID_SITES)
        n = generate_adsorbed(metals_dir, vasp_dir, args.slab, sites, args.calculation,
                              ldau, kpoints, args.encut, args.sigma, freeze_z, args.overwrite)
    print(f"generated {n} VASP input sets")


if __name__ == "__main__":
    main()
