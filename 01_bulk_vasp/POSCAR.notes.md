# POSCAR explanation

Only line 1 is used for a comment. Additional explanation stays in this file
to keep the requested POSCAR format clean. [VASP POSCAR](https://vasp.at/wiki/POSCAR)
defines the scale, lattice rows, species/counts and coordinate block.

| POSCAR line | Meaning and source |
| --- | --- |
| 1 | Structure identification; 12 atoms of LiCoO2 bulk. |
| 2 | Scale 1.0: QE CELL_PARAMETERS are already in angstrom. |
| 3-5 | Three lattice vectors copied token-for-token from QE; no rotation, rescaling or symmetrization. |
| 6 | Li Co O, exactly the QE ATOMIC_SPECIES order and the required POTCAR order. |
| 7 | 3 3 6, counted from the twelve QE atom rows. |
| 8 | Direct: QE ATOMIC_POSITIONS (crystal) are fractional coordinates. |
| 9-11 | Li sites in original order. |
| 12-14 | Co sites in original order. |
| 15-20 | O sites in original order. |

Every numeric token is retained, including negative zeros and the second Co
site (POSCAR line 13) at `1.0000000000 1.0000000000 0.5000000000`.
That site is periodically equivalent to (0,0,0.5), but it has deliberately
not been wrapped. Fractional coordinates on the cell boundary are valid.
No extra atom is created by this boundary representation.

There are no selective-dynamics restrictions in the QE input, so none are
introduced. All ions can relax; INCAR ISIF=3 also enables cell shape/volume.
The file is an initial geometry, not a relaxed CONTCAR. Atomic masses are
carried into INCAR POMASS, since POSCAR contains no mass fields.
