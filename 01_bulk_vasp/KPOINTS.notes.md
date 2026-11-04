# KPOINTS explanation

Only the first line contains a descriptive comment. The remaining four
lines are kept as plain input records. See [VASP KPOINTS](https://vasp.at/wiki/KPOINTS).

| Line | Value | Reason |
| --- | --- | --- |
| 1 | QE 12x12x4 unshifted - Gamma mesh | Describes the original grid and centering. |
| 2 | 0 | Automatic generation from the following grid dimensions. |
| 3 | Gamma | Includes Gamma and corresponds to the unshifted QE mesh. |
| 4 | 12 12 4 | Same subdivisions as QE K_POINTS (automatic). |
| 5 | 0 0 0 | No additional fractional mesh shift. |

QE specifies `12 12 4 0 0 0`: its last three integers are zero shift flags.
The full fractional reciprocal grid consists of (i/12,j/12,k/4), modulo
reciprocal lattice translations, for i,j=0..11 and k=0..3: 576 points.
Symmetry reduction can produce a smaller irreducible set; that number has
not been predicted or assumed to match QE's reduction.

For even subdivisions a VASP Monkhorst-Pack grid would have a different
offset. Gamma centering preserves the supplied grid. This is a mesh with
many k points, not a Gamma-only calculation: use vasp_std, not vasp_gam.
KPAR=4 divides MPI ranks into k-point groups; it does not change this mesh.
The grid is retained as requested and has not been convergence-tested.
