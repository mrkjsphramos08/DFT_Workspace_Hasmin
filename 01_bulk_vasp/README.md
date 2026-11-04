# Annotated LiCoO2 bulk VASP inputs

This package converts the supplied 12-atom `LCO_bulk.in` into a PBE+U+D3
full-cell relaxation setup. It retains every lattice and fractional-position
numeric token from the source, including Co at (1,1,0.5). It uses the requested
550 eV cutoff, Gamma-centered 12x12x4 mesh, Co Ueff=3.2 eV, SIGMA=0.13606 eV,
EDIFFG=-0.0103 eV/angstrom, and KPAR=4 for 40 MPI ranks.

This is an annotated starting setup, not a demonstrated converged calculation.
The main method differences are PAW versus USPP, Hubbard projectors, magnetic
seeds, relaxation/stopping rules, and the D3 three-body term described below.

## Files

| File | Purpose |
| --- | --- |
| INCAR | Every explicit setting has VASP-safe # comments tied to QE or identified as a VASP choice. |
| POSCAR, POSCAR.notes.md | Exact supplied geometry, with explanations outside the input records. |
| KPOINTS, KPOINTS.notes.md | Unshifted Gamma-centered mesh and line-by-line explanation. |
| run_vasp.sh | Annotated one-node, 40-MPI-rank Intel oneAPI SLURM script. |
| prepare_potcar.sh | Optional user-run commands using only the confirmed WSL library paths. |
| reference/LCO_bulk.in | Byte-for-byte copy of the supplied QE input for provenance. |
| verify_inputs.py | Read-only standard-library validator; run with Python 3. |
| VERIFICATION.md | Checks actually performed and their limits. |
| MANIFEST.sha256 | Checksums for package files except the manifest itself. |

No POTCAR, UPF contents, VASP executable or license is included. The original
file in Downloads was read only and was not changed. Generated text inputs and
scripts use ASCII with Unix LF endings and no BOM; the reference copy retains
its original bytes. POSCAR and KPOINTS have only a first-line header comment.

## Precise QE to VASP mapping

Conversions use 1 Ry = 13.605693122994 eV and 1 bohr = 0.529177210903 angstrom.
Retained numerical values do not imply identical basis sets or algorithms.

| QE source | VASP setting or treatment | Status |
| --- | --- | --- |
| calculation='vc-relax'; cell_dynamics='bfgs' | ISIF=3, IBRION=2 | Ions, shape and volume can change; VASP CG replaces QE BFGS. |
| ion_dynamics='bfgs' | IBRION=2 | Optimization intent retained; trajectory and step count differ. |
| nstep=100 | NSW=100 | Same maximum ionic-step count. |
| ibrav=0; CELL_PARAMETERS (angstrom) | POSCAR scale 1.0, same three rows | Exact token preservation. |
| nat=12; ntyp=3; ATOMIC_SPECIES | Li Co O, counts 3 3 6 | Same species and atom order. |
| ATOMIC_POSITIONS (crystal) | Direct coordinates | Exact token preservation, with no wrapping. |
| Li 6.941; Co 58.933; O 15.999 | POMASS=6.941 58.933 15.999 | Same per-species masses; not an MD run. |
| li_pbe_v1.4.uspp.F.UPF; co_pbe_v1.2.uspp.F.UPF; o_pbe_v1.2.uspp.F.UPF | GGA=PE; planned PBE64 Li_sv, Co_pv, O | Requested PBE PAW replacement; UPF contents and valence configurations were not inspected. |
| ecutwfc=40.4 Ry | ENCUT=550 eV | 549.670002169 eV, rounded upward by about 0.330 eV. |
| ecutrho=500 Ry | No direct VASP tag | 6802.846561497 eV is a unit conversion only; do not set ENCUT or ENAUG to that value as a mapping. |
| occupations='smearing'; smearing='gaussian' | ISMEAR=0 | Same Gaussian scheme. |
| degauss=0.01 Ry | SIGMA=0.13606 eV | 0.136056931230 eV before rounding. |
| nspin=2 | ISPIN=2 | Collinear spin polarization. |
| starting_magnetization(2)=0.1 | MAGMOM=3*0.0 3*1.0 6*0.0 | Positive Co seeds of 1 muB chosen; not a numerical conversion of 0.1. |
| HUBBARD (ortho-atomic); U Co-3d 3.2 | LDAU=T, LDAUTYPE=2, LDAUL=-1 2 -1, LDAUU=0 3.2 0, LDAUJ=0 0 0 | Co d-shell Dudarev Ueff=U-J=3.2 eV; different projector definition. |
| vdw_corr='grimme-d3' | IVDW=11 | Zero-damping D3 analogue; pairwise implementation differs from QE default three-body treatment. |
| forc_conv_thr=4D-4 Ry/bohr | EDIFFG=-0.0103 eV/angstrom | 0.010284413495 before rounding; negative VASP sign selects forces, not negative physical forces. |
| electron_maxstep=500 | NELM=500 | Same maximum electronic iterations per ionic step. |
| mixing_beta=0.2; mixing_mode='local-TF' | VASP default mixer | No exact mapping; AMIX=0.2 is not asserted equivalent. |
| diagonalization omitted | ALGO=Normal | Blocked Davidson is a VASP-specific analogue. |
| conv_thr omitted | EDIFF=1E-6 eV | Deliberately chosen SCF tolerance, not a unit conversion; see below. |
| K_POINTS automatic: 12 12 4 0 0 0 | Gamma; 12 12 4; 0 0 0 | Same full unshifted reciprocal grid. |
| press omitted | PSTRESS=0 kbar | Explicit zero-pressure target. |
| cell_dofree omitted | ISIF=3; no selective dynamics | All cell degrees allowed, subject to detected symmetry. |
| pseudo_dir='/home/mark/pseudo_library/' | Manual POTCAR preparation | QE's path does not select VASP potentials. |
| outdir='./tmp'; prefix='lcobulk' | Run in SLURM_SUBMIT_DIR; SYSTEM label | VASP uses fixed filenames rather than QE's prefix/output layout. |
| Fresh QE input, no restart supplied | ISTART=0; ICHARG=2 | Start VASP from atomic charge; no cross-code restart conversion. |
| No corresponding QE switch | PREC=Accurate; LREAL=F; LASPH=T; LMAXMIX=4 | VASP PAW accuracy/mixing choices for this small d-electron cell. |
| User requests 40 MPI ranks | KPAR=4; NCORE=1 | Four groups of ten ranks, with band parallelism inside each group. |

The supplied file does not explicitly set `input_dft`; PBE is inferred from
the source pseudopotential filenames and is also an explicit user request.
Only those names, not their actual UPF headers, have been verified.

## Scientific differences to retain in any report

**Hubbard correction.** The number 3.2 eV is retained, but QE ortho-atomic
orbitals and VASP PAW on-site projectors define different occupation matrices.
The same U value need not yield identical physics across the two codes.
VASP [LDAUTYPE=2](https://vasp.at/wiki/LDAUTYPE) uses U-J; J=0 is the choice
made here. LMAXMIX=4 retains the d-shell density contributions in PAW mixing.

**Dispersion.** With the omitted QE D3 options, documented defaults are
zero damping and three-body enabled. VASP [IVDW=11](https://vasp.at/wiki/DFT-D3)
uses pairwise D3(zero), so this package is not a complete numerical replica.
The source QE version/output is unavailable; the defaults refer to the
[QE 7.5 input documentation](https://www.quantum-espresso.org/Doc/INPUT_PW.html).
A future comparison could consistently omit the three-body term in QE, or
use a validated three-body-capable VASP build. Neither change was made here.
VASP [simple-DFT-D3](https://vasp.at/wiki/Simple-DFT-D3) supports ATM terms
through VDW_S9, requires VASP >=6.6.0 and a build linked to that library, and
can also differ in cutoff treatment. It is not enabled in this INCAR.

**Magnetic seed.** QE 0.1 means polarization per valence electron, whereas
VASP [MAGMOM](https://vasp.at/wiki/MAGMOM) specifies muB per atom. The source
Co UPF valence count was not checked. A +1 muB Co seed retains a positive
spin-polarized starting intent; it neither fixes moments nor establishes
the ground state. Test relevant magnetic initializations before comparing
energies; no constrained total spin or spin-orbit coupling is introduced.

**Convergence and optimization.** QE's omitted conv_thr defaults to 1e-6 Ry
as an estimated SCF error; EDIFF=1e-6 eV instead measures VASP energy changes.
QE tests force components, while VASP's negative EDIFFG tests force norms.
QE also has separate ionic energy and pressure tests (omitted defaults
1e-4 Ry and 0.5 kbar). They are not reproduced by a single EDIFFG value.
Inspect residual stress independently; the requested rounded force bound
does not certify pressure convergence. Compare [EDIFF](https://vasp.at/wiki/EDIFF)
and [EDIFFG](https://vasp.at/wiki/EDIFFG). Hitting NSW or NELM is not success.

**Cell stress and cutoff.** ISIF=3 allows full cell relaxation, but default
symmetry can retain the initial space-group constraints. Replacing USPP
with PAW and keeping about the same cutoff does not establish convergence.
VASP's [ISIF guidance](https://vasp.at/wiki/ISIF) calls attention to Pulay
stress and suggests testing about 1.3 times the largest dataset ENMAX for
cell-volume changes. Keep the requested 550 eV here, then test higher
cutoffs and the k grid for forces, stress and lattice parameters. Local
ENMAX values were not read; 550 eV is not certified sufficient for them.

## POTCAR: optional manual preparation in Ubuntu WSL

Existence of these nonempty files was confirmed in the user's Ubuntu-26.04
distribution on 2026-10-28. Their contents were not opened or changed:

```text
/home/mjram/POTCAR/PBE64/Li_sv/POTCAR
/home/mjram/POTCAR/PBE64/Co_pv/POTCAR
/home/mjram/POTCAR/PBE64/O/POTCAR
```

Li_sv includes Li 1s valence and Co_pv includes Co 3p semicore valence; O is
the standard oxygen choice. These are the requested dataset variants.
See [VASP dataset selection](https://vasp.at/wiki/Choosing_pseudopotentials).
Forty CPU ranks do not determine the functional or potential release.
Dataset header identity, provenance and compatibility still need local review.

From this extracted package directory in WSL, run the supplied commands only
when you are ready to create your licensed local combined file:

```bash
bash prepare_potcar.sh
```

The script rechecks all three files, concatenates Li_sv then Co_pv then O
without inserting comments, and refuses to overwrite POTCAR. It prints
TITEL and ENMAX for your inspection. Confirm the three headers, their order,
the actual release and their cutoff requirements. Codex did not execute it.
Do not add annotations to or distribute the resulting licensed POTCAR.
Keep dataset choices consistent across energies being compared.

The WSL path is local to this computer, not an assumed path on Hasmin.
Transfer the privately prepared POTCAR with the other inputs using your
authorized cluster workflow; no upload, login or submission occurred here.

## Cluster run: 40 MPI ranks, Intel oneAPI, KPAR=4

The script assumes Hasmin partition `high`, one node with 40 available CPUs,
and `/opt/intel/oneapi/setvars.sh`, based on the prior chat. These site
assumptions have not been checked on the cluster. Add its required time,
memory, account or QOS directives before the first executable line.
The script requires a Slurm allocation and refuses an incorrect rank count.

After sourcing oneAPI, `vasp_std` must be on PATH, or set VASP_EXE to its
absolute path before submission. oneAPI initialization does not install or
license VASP. The binary must match the Intel MPI runtime used by mpiexec.
ONEAPI_SETVARS can select a different trusted site environment script.

The launch is `mpiexec -np "$SLURM_NTASKS" "$vasp_exe"`, with one OMP and
MKL thread per rank and Intel Hydra Slurm bootstrap. See
[Intel scheduler support](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-16/job-schedulers-support.html).
Confirm the site's launcher/module policy; no PMI library path is guessed.
Do not run this script directly on a login node.

[KPAR=4](https://vasp.at/wiki/KPAR) divides forty ranks into groups of ten.
NCORE=1 is a starting layout, not a benchmark result. Actual speed, memory,
band rounding and load balance depend on the build and irreducible k points.
Use `vasp_std`: Gamma centering does not make this a Gamma-only job.

After manual POTCAR preparation and site configuration, from the calculation
directory on the cluster:

```bash
python3 verify_inputs.py
bash -n run_vasp.sh
grep TITEL POTCAR
grep ENMAX POTCAR
sbatch run_vasp.sh
```

The validator checks the distributable package and deliberately does not
read POTCAR. After running VASP, check interpreted tags and dataset order
in OUTCAR, SCF convergence at each step, final forces and stress, and the
final structure. A zero launcher exit code is not physical validation.

## Reproducing the package checks

```bash
python3 verify_inputs.py
sha256sum -c MANIFEST.sha256
bash -n run_vasp.sh
bash -n prepare_potcar.sh
```

VERIFICATION.md records what was checked locally. No electronic calculation,
cutoff/k-point convergence study, final magnetic-state validation or cluster
runtime test is claimed.
