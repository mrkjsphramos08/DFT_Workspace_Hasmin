# Verification receipt

Prepared 2026-10-28 from the supplied `LCO_bulk.in`.

## Passed checks

- Original Downloads file and bundled reference are byte-for-byte identical.
- All 9 lattice and 36 fractional-position numeric tokens match the QE source
  exactly, including signed zeros and Co at (1,1,0.5).
- Twelve atoms, Li3 Co3 O6, with original atom/species order and no constraints.
- Cell is right-handed and nonsingular; volume 96.681518153582
  angstrom^3. Independent minimum periodic separation is
  1.928964207465 angstrom.
- POSCAR and KPOINTS contain clean syntax records, with comments only in the
  first-line headers. INCAR has a comment explaining every explicit setting.
- Standard-library validator checks all 31 explicit INCAR values, the Hubbard
  array ordering and Ueff, twelve expanded MAGMOM seeds, zero-shift Gamma
  12x12x4 mesh, Slurm directive ordering and the 40/4=10 rank grouping.
- Independent ASE 3.29.0 readers parsed both QE and POSCAR;
  species, cells and unwrapped fractional positions agree with absolute
  tolerance 1e-13. ASE's VASP readers independently parsed INCAR and KPOINTS,
  checking the principal physics/parallel tags, smearing, force threshold,
  cutoff, magnetization array length, mesh and Gamma flag.
- Both shell scripts pass `bash -n` in Ubuntu-26.04 WSL.
- The run script exits with status 2 for no Slurm allocation, wrong MPI rank
  count, and missing POTCAR. These guard checks stop before sourcing oneAPI
  or invoking MPI. Slurm environment values were simulated locally only.
- Generated text uses LF endings without a BOM. The source reference retains
  its original bytes. Shell files have executable permission metadata in ZIP.
- Archive CRC, member bytes, checksums and validator are checked again after
  extraction into a local scratch directory.

## Scope and limitations

PBE64 Li_sv, Co_pv and O files were confirmed nonempty at the requested
literal WSL paths earlier in this chat. Dataset contents, header identities,
actual ENMAX and licensed release provenance were not inspected. POTCAR was
not read, generated, edited, included or uploaded. prepare_potcar.sh was
syntax-checked only; its concatenation commands were not executed.

No VASP run, MPI launch, remote login, job submission or cluster policy check
was performed. Intel oneAPI path, executable/runtime ABI, memory, wall time,
partition suitability, performance and final scientific convergence remain
unverified. This receipt establishes conversion consistency and offline input
syntax, not successful VASP execution. README.md documents the D3 three-body,
PAW/USPP, Hubbard-projector, magnetic-seed and stopping-rule differences.
