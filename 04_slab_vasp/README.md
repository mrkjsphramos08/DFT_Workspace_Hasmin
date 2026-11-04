# LiCoO2 (104) 3x2 Slab Relaxation (VASP)

## Calculation Overview
- **Structure:** LiCoO2 (104) surface slab ($3 \times 2$ in-plane supercell)
- **Cell Dimensions:** $a = 8.431\ \text{Å}$, $b = 11.356\ \text{Å}$, $c = 22.142\ \text{Å}$ ($\alpha = \beta = 90^\circ, \gamma = 90^\circ$)
- **Composition:** 72 atoms (18 Li, 18 Co, 36 O)
- **Vacuum:** $18.0\ \text{Å}$ placed above the top surface ($z = 4.14 \rightarrow 22.14\ \text{Å}$)
- **Selective Dynamics:**
  - **Bottom 2 layers (48 atoms):** Fixed (`F F F`)
  - **Top layer (24 atoms):** Free to relax (`T T T`)

## Input Files
- `POSCAR`: Generated with selective dynamics from the relaxed bulk structure (`01_vc_relax_vasp_gpt/CONTCAR`).
- `INCAR`: Adapted from bulk PBE+U(3.2)+D3 calculation with key slab modifications:
  - `ISIF = 2`: Relaxes ion positions while keeping cell shape and volume fixed (prevents collapse of vacuum).
  - `MAGMOM = 18*0.0 18*1.0 36*0.0`: $1.0\ \mu_B$ initial seed per Co atom.
  - `LDIPOL = .TRUE.`, `IDIPOL = 3`, `DIPOL = 0.5 0.5 0.0935`: Electrostatic dipole correction along the slab normal ($z$).
  - `KPAR = 2`, `NCORE = 4`: Optimized parallel distribution over k-points and bands on 40 MPI ranks (eliminates single-core FFT bottleneck).
- `KPOINTS`: Gamma-centered $4 \times 3 \times 1$ grid.
- `POTCAR`: PAW PBE datasets in `Li_sv`, `Co_pv`, `O` order (matching POSCAR).
- `run_vasp.sh`: Slurm job submission script (1 node, 40 MPI tasks, `high` partition).
- `python_vasp_analyze.py`: Analysis script to summarize convergence, residual forces, and energies.

## Submission
```bash
cd /home/mark/LCO/01_bulk/04_slab_vasp
sbatch run_vasp.sh
```

