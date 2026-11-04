# `05_metals_vasp/` — VASP mirrors of the Quantum ESPRESSO metal‑cluster campaign

This directory holds **isolated‑cluster** and **adsorbed** VASP input decks that
parallel the existing `05_metals/` (QE) work tree.  VASP uses the **same relaxed
structures** (CIFs) as the QE side, so the two engines can be compared directly.

## Convention

- **Structure provenance**: Every VASP `POSCAR` is generated from a *real* relaxed
  structure saved under `.../<Metal>/.../<cif>`.  The source CIF is kept under the
  VASP tree (`05_metals_vasp/<...>/.../<cif>`) for provenance and round‑tripping.
- **Engine dimension**: `05_metals/` → QE; `05_metals_vasp/` → VASP.  Adding a new
  DFT engine never touches the other’s directory, which keeps git history
  clean and avoids accidental breakage.
- **Files tracked**: CIF (structure) + `INCAR` (parameters) + `KPOINTS` (k‑mesh)
  + `POSCAR` (input).  Everything else (`POTCAR*`, `OUTCAR*`, `vasprun.xml`,
  `OSZIC*`, `WAVECAR`, `CHGCAR`, `tmp/`… ) follows the same git‑ignore policy as
  the QE side.
- **Parameter mapping**: see the root `README.md` for the QE ↔ VASP table.

## Seeds

Only a **minimal subset of systems** is committed to keep the repository size
reasonable and avoid liability for proprietary VASP PAW datasets.  The seeds are
kept deliberately small but complete (ready to run after a user sets
`$VASP_PP_PATH`).

### `01_vc_relax_vasp/`

* Bulk LiCoO₂ (`LCO_bulk.vc_relax.cif`) – PBE+U(Co)=3.2 eV + D3 zero‑damping,
  vc‑relax (IBRION=1, ISIF=3).  KPOINTS: 12×12×4 MP.

### `05_metals_vasp/Co4/Tetra/`

* Bare Co₄ tetrahedron (`Co4_Tetra_latest.cif`) – PBE+D3, Γ‑point, pure
  magnetic moment of 0.80 μB per atom (matches QE `starting_magnetization=0.8`).
  Relax (`IBRION=2, ISIF=2, NSW=200`).  No Hubbard U.

> **All other metals** (Cu₄, Ni₄, Pt₄, Pt₁) **can be generated** by the
> `vasp/generate_inputs.py` script described below.

## Using the VASP generators

The whole campaign (clusters + adsorbed systems) can be reproduced in VASP
with a single Python call.  No hand‑editing is required after the initial
seed — all INCARs, KPOINTs and POSCARs are produced from the existing relaxed
CIFs.

### Prerequisites

```bash
# Minimal test setup (without real VASP)
python -m vasp.generate_inputs --all

# Only a specific metal, e.g. Pt4 clusters & adsorbed Pt1
python -m vasp.generate_inputs --metal Pt4 --cif <path>
```

### CLI synopsis

```bash
python -m vasp.generate_inputs \
  --cif <path/to/relaxed.cif> \
  --out-dir <work/vasp> \
  [--calculation {vc-relax,relax,scf}] \
  [--kpoints 4 3 1] \
  [--ldaU "Co=3.2,Ni=6.0"] \
  [--freeze-below 0.5] \
  [--overwrite]
```

### Key options

* `--calculation` — `vc-relax`, `relax` (default), or `scf`.
* `--kpoints` — MP grid (default 4×3×1 for surfaces, 12×12×4 for bulk).
* `--ldaU` — space‑separated list `Element=U` (use `"off"` to skip U).
* `--freeze-below` — freeze bottom ``z`` of a slab (affects adsorbed systems).
* `--no-freeze` — turn off automatic slab freezing.

#### Metal‑campaign shortcuts

```bash
# All bare clusters (no U)
python -m vasp.generate_inputs --clusters

# All adsorbed systems (Co, Ni, Cu, Pt) – frozen slab, Co U=3.2
python -m vasp.generate_inputs --slab 3_layers --site li,co,lico,center,oxygen,lili,coco,lio,coo

# One specific system, e.g. Co4/Tetra/Absorbed/3_layers/li
python -m vasp.generate_inputs --metal Co4 --shape Tetra --slab 3_layers --site li
```

The script writes:

```
<out-dir>/INCAR
<out-dir>/KPOINTS
<out-dir>/POSCAR
<out-dir>/<cif>
<out-dir>/POTCAR_SPECIES
```

and prints the relative path for easy navigation.

## SLURM scripts

`vasp/run_sh.py` mirrors the QE `run_sh.*` tools:

```bash
python -m vasp.run_sh --vasp-dir 05_metals_vasp
```

It creates a ``run_vasp.sh`` in every directory that now contains an
``INCAR`` but no ``OSZIC``/`OUTCAR`, so each experiment is re‑run only
once.  Each script:

* uses ``module load vasp/<version>`` (user‑configurable),
* assembles ``POTCAR`` from ``$VASP_PP_PATH`` using the ``POTCAR_SPECIES`` file,
* runs ``mpirun -n $SLURM_NTASKS vasp_std``.

## VASP POTCAR handling

PAW datasets are **not** stored in the repository.  The generated
``POTCAR_SPECIES`` lists the element‑order needed to concatenate the user’s
PBE PAW library:

```bash
cd 05_metals_vasp/Co4/Tetra/Adsorbed/3_layers/lico
python ../../../vasp/setup_potcar.py \
  --pp-path /my/pbe/potcar \
  --subdir PBE \
  --species-file POTCAR_SPECIES
```

This produces ``POTCAR`` in the current work directory.

## Validation and CI

* `lco/vasp_output.py` parses ``OUTCAR``/``OSZIC`` (pure‑text, no VASP required).
* `tests/test_vasp_output.py` contains unit tests that run in CI (no ASE needed).
* The CI also imports the VASP modules and runs `python -m compileall .`

## Getting started

1. Add your VASP PBE PAW library path (`$VASP_PP_PATH`).
2. Choose an engine: **parallel** `05_metals_vasp/` (recommended) or **inline**
   `<metal>/Tetra/Adsorbed/<slab>/<site>/` (if you are comfortable with
   mixing QE and VASP in the same directories).
3. Use `vasp/generate_inputs.py` to create all VASP input decks.
4. Run `vasp/run_sh.py` to emit ``run_vasp.sh`` scripts.
5. Assemble ``POTCAR`` with `vasp/setup_potcar.py` and submit with ``sbatch``.

---

**Quick cheat‑sheet** (edit `dft_terminal_cheatsheet.md` to add VASP lines):

```bash
# generate all VASP inputs from QE‑relaxed structures
python -m vasp.generate_inputs --all

# assemble POTCAR (one per site, script uses POTCAR_SPECIES)
python -m vasp.setup_potcar.py --pp-path <path>

# submit one job
cd some/vasp/dir
bash run_vasp.sh
```

Enjoy the parallel VASP campaign!