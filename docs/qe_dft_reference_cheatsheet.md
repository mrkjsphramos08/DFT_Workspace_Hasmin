# Quantum ESPRESSO DFT Reference Cheat Sheet

A non-terminal reference for building and understanding Quantum ESPRESSO (QE) `pw.x` input
files — namelists, cards, and commonly-used keywords, with notes on Hubbard U, vdW corrections,
and spin polarization.

## How to use this file

- This is a **reference for writing/editing `.pwi` input files**, not terminal commands (see the
  companion terminal cheat sheet for that).
- Every QE `pw.x` input file is built from **namelists** (`&CONTROL`, `&SYSTEM`, `&ELECTRONS`,
  optionally `&IONS` and `&CELL`) followed by **cards** (`ATOMIC_SPECIES`, `ATOMIC_POSITIONS`,
  `K_POINTS`, `CELL_PARAMETERS`).
- Each namelist ends with a line containing only `/`.
- Skim the section headers below to find the block you need, copy the relevant keyword lines,
  and edit values to match your system.

---

## 1. Overall Input File Skeleton

```
&CONTROL
  calculation = 'relax'
  ...
/
&SYSTEM
  ibrav = 0
  nat = 20
  ntyp = 3
  ecutwfc = 44.1
  ecutrho = 441
  ...
/
&ELECTRONS
  ...
/
&IONS
  ...
/
&CELL
  ...
/
ATOMIC_SPECIES
  ...
ATOMIC_POSITIONS {angstrom}
  ...
K_POINTS {automatic}
  ...
CELL_PARAMETERS {angstrom}
  ...
```

`&IONS` is only needed for `calculation = 'relax'/'md'/'vc-relax'`; `&CELL` only for
`'vc-relax'/'vc-md'`.

---

## 2. `&CONTROL` Namelist — Job Setup

| Keyword | Purpose | Common values |
|---|---|---|
| `calculation` | Type of run | `'scf'`, `'relax'`, `'vc-relax'`, `'md'`, `'nscf'`, `'bands'` |
| `restart_mode` | Fresh start or restart | `'from_scratch'`, `'restart'` |
| `prefix` | Label for output files | any string, e.g. `'lco_dopant'` |
| `outdir` | Scratch/output directory | e.g. `'./tmp/'` |
| `pseudo_dir` | Path to pseudopotential files | e.g. `'./pseudo/'` |
| `tprnfor` | Print forces on atoms | `.true.` |
| `tstress` | Print stress tensor | `.true.` |
| `forc_conv_thr` | Force convergence threshold (Ry/au) | `1.0e-3` |
| `etot_conv_thr` | Energy convergence threshold (Ry) | `1.0e-4` |
| `nstep` | Max ionic/electronic steps | e.g. `100` |

---

## 3. `&SYSTEM` Namelist — Structure & Physics

| Keyword | Purpose | Common values |
|---|---|---|
| `ibrav` | Bravais lattice index | `0` (use `CELL_PARAMETERS` explicitly) |
| `nat` | Number of atoms | integer |
| `ntyp` | Number of atomic species | integer |
| `ecutwfc` | Wavefunction cutoff (Ry) | e.g. `44.1` |
| `ecutrho` | Charge density cutoff (Ry) | typically `4×–10×` ecutwfc for ultrasoft/PAW; e.g. `441` |
| `occupations` | Occupation scheme | `'smearing'`, `'fixed'` |
| `smearing` | Smearing type | `'gaussian'`, `'mp'`, `'mv'` |
| `degauss` | Smearing width (Ry) | e.g. `0.01` |
| `nspin` | Spin polarization | `1` (none), `2` (collinear) |
| `starting_magnetization(i)` | Initial magnetic moment per species `i` | e.g. `0.5` |
| `nosym` | Disable symmetry (useful for doped/defect cells) | `.true.` |

### Hubbard U (DFT+U) — modern QE syntax (v6.7+)

```
&SYSTEM
  ...
  lda_plus_u = .true.
  lda_plus_u_kind = 0
/
HUBBARD {ortho-atomic}
  U Co-3d  4.0
```

- Older QE versions instead use `Hubbard_U(i)` inside `&SYSTEM`, matched to species index `i`.
- `lda_plus_u_kind = 0` is simplified (Dudarev) DFT+U; `= 1` is full rotationally-invariant.

### van der Waals correction (DFT-D3)

```
&SYSTEM
  ...
  vdw_corr = 'dft-d3'
  dftd3_version = 3
```

- `dftd3_version = 2` (Grimme D2), `3` (D3 zero-damping), `4` (D3 Becke-Johnson damping).

---

## 4. `&ELECTRONS` Namelist — SCF Convergence

| Keyword | Purpose | Common values |
|---|---|---|
| `conv_thr` | SCF energy convergence threshold (Ry) | e.g. `1.0e-8` |
| `mixing_mode` | Charge mixing scheme | `'plain'`, `'local-TF'` |
| `mixing_beta` | Mixing factor | `0.3`–`0.7` (lower = more stable, slower) |
| `electron_maxstep` | Max SCF iterations | e.g. `100` |
| `diagonalization` | Diagonalization algorithm | `'david'`, `'cg'` |

---

## 5. `&IONS` Namelist — Ionic Relaxation

| Keyword | Purpose | Common values |
|---|---|---|
| `ion_dynamics` | Relaxation algorithm | `'bfgs'`, `'damp'` |
| `upscale` | BFGS trust radius scaling | e.g. `100` |

---

## 6. `&CELL` Namelist — Variable-Cell Relaxation

| Keyword | Purpose | Common values |
|---|---|---|
| `cell_dynamics` | Cell relaxation algorithm | `'bfgs'` |
| `press` | Target external pressure (kbar) | `0.0` for ambient |
| `cell_dofree` | Constrain which cell DOFs relax | `'all'`, `'2Dxy'` (useful for slabs) |

---

## 7. Cards

### `ATOMIC_SPECIES`
```
ATOMIC_SPECIES
Li  6.941   Li.pbe-sl-rrkjus_psl.1.0.0.UPF
Co  58.933  Co.pbe-spn-rrkjus_psl.1.0.0.UPF
O   15.999  O.pbe-n-rrkjus_psl.1.0.0.UPF
```
Format: `symbol  atomic_mass  pseudopotential_filename`

### `ATOMIC_POSITIONS`
```
ATOMIC_POSITIONS {angstrom}
Co   0.000  0.000  0.000
O    1.418  1.418  1.418  0  0  1
```
- Units: `{angstrom}`, `{crystal}` (fractional), `{bohr}`.
- Optional trailing `0 1 1` flags fix (`0`) or free (`1`) an atom along x/y/z — useful for fixing
  slab bottom layers during relaxation.

### `K_POINTS`
```
K_POINTS {automatic}
4 3 1  0 0 0
```
- Format: `nk1 nk2 nk3  shift1 shift2 shift3` (shift = `0` or `1`, i.e. Monkhorst-Pack shifted or not).
- For slabs, keep the out-of-plane (vacuum-direction) k-point at `1`.
- Alternative modes: `{gamma}` (Γ-point only, no shift line needed), `{crystal_b}` (band structure path).

### `CELL_PARAMETERS`
```
CELL_PARAMETERS {angstrom}
  a1x  a1y  a1z
  a2x  a2y  a2z
  a3x  a3y  a3z
```
Units: `{angstrom}`, `{bohr}`, or `{alat}`.

---

## 8. Common Convergence / Troubleshooting Cheatsheet

| Symptom | Likely fix |
|---|---|
| SCF not converging | Lower `mixing_beta` (e.g. to `0.1`–`0.3`); increase `electron_maxstep` |
| Ionic relaxation stalls | Loosen `forc_conv_thr`; check `upscale`/BFGS trust radius |
| Imaginary/negative charge density warnings | Increase `ecutrho`; check pseudopotential compatibility |
| Metallic system not converging | Use `occupations='smearing'` with appropriate `degauss` |
| Doped/defect cell symmetry errors | Set `nosym = .true.` |
| Slab drifting/rotating during `vc-relax` | Use `cell_dofree = '2Dxy'` to keep vacuum direction fixed |

---

### Notes

- Namelist keywords are case-insensitive; card names (`ATOMIC_SPECIES`, etc.) are conventionally
  uppercase.
- `ecutrho` should always be checked against your pseudopotential type — norm-conserving PPs
  often need only `4×ecutwfc`, while ultrasoft/PAW can need `8×–12×`.
- Full authoritative keyword list: the official [PW.x input description](https://www.quantum-espresso.org/Doc/INPUT_PW.html) is the ground truth if a keyword's default behavior is unclear.
