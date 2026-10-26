# LCO Metal-Cluster Adsorption — DFT Workflow

Automated [Quantum ESPRESSO](https://www.quantum-espresso.org/) workflow for studying
**M₄ metal clusters adsorbed on LiCoO₂ (LCO) surfaces**, from bulk relaxation through
slab construction, cluster generation, site adsorption, job submission and analysis.

[![CI](https://github.com/mrkjsphramos08/DFT_Workspace_Hasmin/actions/workflows/ci.yml/badge.svg)](https://github.com/mrkjsphramos08/DFT_Workspace_Hasmin/actions/workflows/ci.yml)

Everything is driven by small, composable Python scripts built on
[ASE](https://wiki.fysik.dtu.dk/ase/), with shared logic factored into the in-repo
`lco/` package, so each stage can be re-run independently and the whole campaign is
reproducible from the inputs kept in this repository.

| | |
|---|---|
| Code | [Quantum ESPRESSO](https://www.quantum-espresso.org/) **7.2** (`pw.x`), MPI on SLURM |
| Functional | PBE + Hubbard U (`U Co-3d 3.2 eV`, ortho-atomic) + Grimme **D3** dispersion |
| Spin | `nspin = 2`, gaussian smearing (`degauss = 0.01 Ry`) |
| Pseudopotentials | PBE ultrasoft (`.uspp.F.UPF`): Li, Co, O, Pt, Ni, Cu |
| Cutoffs | `ecutwfc = 40.4 Ry`, `ecutrho = 500 Ry` |
| Slab | 3×2 LiCoO₂ surface, 3-layer (72 atoms) and 4-layer (96 atoms), ~24 Å vacuum |
| Clusters | M₄ tetrahedra (Pt₄, Co₄, Cu₄, Ni₄) and rhombi in a 15 Å cubic box |
| Python | 3.10 · ASE 3.29 · NumPy 2.2 · pandas 2.3 |

---

## Repository layout

```
.
├── 01_vc_relax/            bulk LiCoO₂ variable-cell relaxation
├── 02_scf/                 placeholder — fixed-cell SCF follow-up
├── 03_dos/                 placeholder — DOS analysis
├── 04_slab/
│   ├── 3_layers/           3×2 LCO slab, 3 layers  (72 atoms)
│   ├── 4_layers/           3×2 LCO slab, 4 layers  (96 atoms)
│   ├── LCO_slab_distance.py    in-plane layer-distance analysis
│   └── lco_slab_distances.csv  exported distances
├── 05_metals/              clusters + adsorbed systems
│   └── <Metal>/            Co4 · Cu4 · Ni4 · Pt4
│       ├── Rhom/           rhombic isomer
│       └── Tetra/          tetrahedral isomer
│           └── Adsorbed/
│               ├── 3_layers/<site>/   slab-specific adsorption set
│               └── 4_layers/<site>/
├── 05_metals_legacy/       archived 44.1/441 generation (NOT tracked in git)
├── Adsorbed_clusters_tetra_new.py
├── Adsorbed_QE_lico.py
├── adsorption_config.py
├── adsorption_sites.py
├── archive_legacy_cluster_runs.py
├── cleanup_adsorbed_folders.py
├── migrate_adsorption_folders.py
├── prepare_cluster_rerun.py
├── regenerate_cif.py
├── run_sh.py
├── run_sh_clusters.py
├── Surface_energy.py
├── dft_terminal_cheatsheet.md
└── qe_dft_reference_cheatsheet.md
```

Raw Quantum ESPRESSO outputs, SLURM logs and `tmp/` scratch are **not tracked** —
they are regenerable and total ~17 GB here. Only code, documentation, QE inputs and
reference structures are committed.

---

## Naming conventions

Understanding the file names is the key to reading this tree.

| pattern | meaning |
|---|---|
| `INPUT_<name>.in` / `OUTPUT_<name>.out` | pw.x input and the output it produced |
| `..._ecut40.4_500.in` | input re-generated at `ecutwfc = 40.4` / `ecutrho = 500` |
| `<Metal>_<Shape>_latest.cif` | structure extracted from the newest relaxed output |
| `adsorbed_<cluster-stem>_<site>.cif` | adsorbed structure at one site |
| `Adsorbed/{3_layers,4_layers}/<site>/` | adsorption namespace per slab profile |
| `05_metals_legacy/` | mirror of the superseded 44.1/441 generation |

The `_ecut40.4_500` tag is what distinguishes the **current** generation from the
archived 44.1/441 one. When the cutoffs changed, `prepare_cluster_rerun.py` wrote
tagged copies rather than overwriting, so both generations coexisted until
`archive_legacy_cluster_runs.py` moved the old ones into `05_metals_legacy/`.

---

## Workflow

The numbered directories are the stages, in order. Each stage consumes the previous
one and is driven by a script that can be re-run at any time.

### 1 — Bulk relaxation (`01_vc_relax`)

Variable-cell relaxation of bulk LiCoO₂, giving the lattice constants used to build
the slabs. The `U = 3.2 eV` Hubbard value is set here and inherited downstream.

```bash
cd 01_vc_relax && sbatch run.sh          # writes LCO_bulk.out
```

### 2 — Slab construction and relaxation (`04_slab`)

A 3×2 surface cell is cleaved from the relaxed bulk at 3- and 4-layer thicknesses and
relaxed with the lower layers frozen. `adsorption_config.py` records which CIF each
profile uses:

| profile | slab file | atoms | output namespace |
|---|---|---|---|
| `3_layers` | `04_slab/3_layers/3x2_3_layers.cif` | 72 | `Adsorbed/3_layers/` |
| `4_layers` | `04_slab/4_layers/3x2_4_layers.cif` | 96 | `Adsorbed/4_layers/` |

Both CIFs are the final relaxed geometries of their own runs and both use the current
40.4/500 cutoffs. `DEFAULT_SLAB = "3_layers"`.

`04_slab/LCO_slab_distance.py` measures adjacent in-plane distances per layer and
writes `lco_slab_distances.csv`.

### 3 — Isolated clusters (`05_metals/<Metal>/{Tetra,Rhom}`)

Each metal is relaxed as an isolated M₄ cluster in a 15 Å box, in two geometries:

| cluster | tetrahedron edge (Å) |
|---|---|
| Co₄ | 2.2548 |
| Cu₄ | 2.3882 |
| Ni₄ | 2.2738 |
| Pt₄ | 2.5835 |

Full bond-length tables are in `metal_cluster_bond_lengths.csv`; energies and
magnetisations in `metal_clusters_specifications.csv`. Bare clusters use Γ-point
sampling.

### 4 — Adsorption (`05_metals/<Metal>/Tetra/Adsorbed/<profile>/<site>/`)

Seven sites are generated per cluster per slab profile: `li`, `co`, `lico`,
`centered`, `oxygen`, `lili`, `coco`. The cluster is oriented base-down and placed
2.0 Å above the topmost surface atoms, and the bottom 50 % of the slab height is
frozen — for a 3-layer slab that is the bottom 2 of 3 layers — via `0 0 0` flags in
`ATOMIC_POSITIONS`.

```bash
# 1. build the adsorbed structures
python Adsorbed_clusters_tetra_new.py --slab 3_layers

# 2. build the pw.x inputs from the cluster template
python Adsorbed_QE_lico.py --slab 3_layers

# 3. create SLURM scripts for sites that have no output yet
python run_sh.py --slab 3_layers
```

### 5 — Submission and analysis

`run.sh` must be submitted **from inside its own site directory**, because the
scripts use relative paths and `outdir = './tmp/'`:

```bash
cd 05_metals/Pt4/Tetra/Adsorbed/3_layers/oxygen
sbatch run.sh
```

When a run finishes, refresh its structure and analyse:

```bash
python regenerate_cif.py        # run from the repository root
python Surface_energy.py        # slab/bulk surface energy
```

`run_sh.py` and `Adsorbed_QE_lico.py` deliberately **skip any site directory that
already contains a `.out` file**, so re-running them is safe and idempotent — but it
also means a finished site must have its output moved aside before it can be
regenerated.

---

## Script reference

| script | purpose |
|---|---|
| `adsorption_config.py` | registry of slab profiles (slab CIF + output namespace); `DEFAULT_SLAB` |
| `adsorption_sites.py` | canonical adsorption-site names and their order |
| `Adsorbed_clusters_tetra_new.py` | place relaxed clusters on the slab sites and write `adsorbed_*_<site>.cif` |
| `Adsorbed_QE_lico.py` | build pw.x inputs for adsorbed systems from the cluster template |
| `run_sh.py` | SLURM scripts for adsorbed runs (16 tasks) |
| `run_sh_clusters.py` | SLURM scripts for bare cluster runs (8 tasks) |
| `prepare_cluster_rerun.py` | regenerate cluster inputs at new cutoffs from a relaxed output |
| `regenerate_cif.py` | rebuild `*_latest.cif` from QE outputs |
| `archive_legacy_cluster_runs.py` | move the superseded 44.1/441 generation into `05_metals_legacy/` |
| `migrate_adsorption_folders.py` | move unlabelled `Adsorbed/<site>/` into the `4_layers` namespace |
| `cleanup_adsorbed_folders.py` | archive orphaned CIFs and stale `run.sh` leftovers |
| `Surface_energy.py` | surface energy from the slab and bulk outputs |
| `04_slab/LCO_slab_distance.py` | in-plane layer distances → `lco_slab_distances.csv` |

Common flags: `--slab {3_layers,4_layers}`, `--apply` (the housekeeping scripts are
dry-run by default), `--cluster` / `--shape` on `prepare_cluster_rerun.py`.

### The `lco/` shared library

Shared logic that was previously copy-pasted across scripts lives in the `lco/` package:

| module | provides |
|---|---|
| `lco.paths` | `PROJECT_ROOT`, `PSEUDO_DIR`, `SCRATCH_DIR` — single source of truth for filesystem locations |
| `lco.qe_input` | `PSEUDO_DB`, `qe_command`, `variant_prefix` — building `pw.x` input decks |
| `lco.qe_output` | `parse_scf_history`, `final_accuracy` — parsing QE outputs for convergence analysis |

Import them in any script:

```python
from lco.paths import PROJECT_ROOT, PSEUDO_DIR
from lco.qe_input import PSEUDO_DB, qe_command
```

---

## Testing

The repository ships with a test suite powered by [pytest](https://docs.pytest.org/).

```bash
# run the tests
python -m pytest tests/

# CI performs the same checks:
#   - import all lco submodules
#   - compile every .py file (compileall)
```

`tests/test_qe_output.py` exercises `parse_scf_history` against a sample QE output
shipped in the repository. If the sample is absent the test is auto-skipped.

---

## Environment

```bash
# conda environment used for development
conda create -n dft_clean python=3.10
conda activate dft_clean
pip install -r requirements.txt
```

`requirements.txt` records the tested versions (ASE 3.29.0, NumPy 2.2.6, pandas
2.3.3). On the cluster, Quantum ESPRESSO is provided by a module:

```bash
module load qe/7.2
```

SLURM resources used per job class:

| job class | partition | tasks | example |
|---|---|---|---|
| bare cluster | `high` | 8 | `Co4_Tetra` |
| adsorbed system | `high` | 16 | `Pt4_Tetra_adsorbed_oxygen` |
| slab | `high` | 32 | `3x2_4_layers` |

---

## Gotchas worth knowing

These are the things that actually cost time, written down so they don't again.

1. **Shared `prefix` between generations.** Old and new inputs of the same system
   often share `prefix` and `outdir = './tmp/'`, so a re-run reuses (and overwrites)
   the previous run's scratch. Keep one generation per directory, or archive the old
   one first — this is exactly what the `_ecut` tagging plus
   `archive_legacy_cluster_runs.py` workflow exists for.

2. **ASE cannot parse some of these QE outputs.** `ase.io.read(..., format='espresso-out')`
   raises `AssertionError` on several truncated/multi-step outputs. On this ASE build
   there is also no `espresso-xml` reader registered, so the XML fallback inside
   `regenerate_cif.py` cannot rescue them. The reliable route is to parse the
   `Begin final coordinates … End final coordinates` block of the `.out` directly.

3. **Slab CIFs wrap across the cell boundary.** The vacuum is inside the cell, so
   some atoms sit at `scaled_z` near 0.997. A naive `max(positions[:, 2])` therefore
   returns a *periodic image*, not the surface. Always work in fractional `z`
   (`Adsorbed_clusters_tetra_new.py` does this via `scaled_z < 0.8`).

4. **`regenerate_cif.py` uses a relative default root** (`"05_metals"`), so it must be
   launched from the repository root — from anywhere else it silently finds nothing.

5. **Template selection is order-dependent.** `Adsorbed_QE_lico.py` picks the first
   non-adsorbed `*.in` it globs in the `Tetra` directory as its template. Keep exactly
   one template there (which is what the `_ecut` tagging plus archiving achieves) or
   the choice is undefined.

---

## Reference results

Isolated-cluster benchmark from the archived 44.1/441 generation
(`metal_clusters_specifications.csv`); the 40.4/500 re-runs reproduce the same
relaxed geometries to within 0.0007 Å.

| system | E_total (Ry) | E_Fermi (eV) | magnetisation (μB/cell) |
|---|---|---|---|
| Co₄ Rhom | −1193.114027 | −3.8215 | 10.00 |
| Co₄ Tetra | −1193.079247 | −3.7410 | 9.88 |
| Cu₄ Rhom | −1614.566157 | −3.8175 | 0.00 |
| Cu₄ Tetra | −1614.501018 | −3.7607 | 2.00 |
| Ni₄ Rhom | −1372.338166 | −3.7195 | 4.00 |
| Ni₄ Tetra | −1372.340964 | −3.2829 | 3.97 |
| Pt₄ Rhom | −841.962028 | −4.6835 | 3.84 |
| Pt₄ Tetra | −841.954527 | −3.6626 | 2.55 |

---

## What is deliberately not in the repository

| excluded | why |
|---|---|
| `**/tmp/`, `*.save/` | QE scratch — regenerable, ~17 GB here |
| `*.out` | QE outputs and SLURM logs — regenerable |
| `*.xyz` | derived structure exports |
| `*.tar.gz` | bulk data snapshots (`all_cifs.tar.gz`, `all_outs.tar.gz`) |
| `05_metals_legacy/` | archived 44.1/441 generation |
| pseudopotentials | live in `pseudo_dir` on the cluster (`/home/mark/pseudo_library/`) |

---

## Docs

Two companion references live in the repository root:

- **`dft_terminal_cheatsheet.md`** — one-line ASE/SLURM/disk commands
- **`qe_dft_reference_cheatsheet.md`** — how to build a `pw.x` input by hand, with
  Hubbard U, vdW and spin-polarisation notes

## License

[MIT](LICENSE) — © 2026 mrkjsphramos08
`archive_legacy_cluster_runs.py` moved the old ones out.