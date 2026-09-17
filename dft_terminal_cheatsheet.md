# DFT / Quantum ESPRESSO Terminal Cheat Sheet

A reference of one-line terminal commands for Quantum ESPRESSO (QE) DFT workflows,
ASE-based structure handling, and SLURM/HPC job management.

## How to use this file

- Everywhere you see **`input`**, replace it with your actual file (e.g. `relax.out`, `LOG.out`, `structure.cif`).
- Everywhere you see **`output`**, replace it with the filename you want to create (e.g. `relaxed.xyz`, `final.vasp`).
- `<jobid>`, `<element>`, `job1/`, `job2/` etc. are placeholders — swap in real values.
- Commands starting with `python -O -c "..."` are single-line Python scripts using [ASE](https://wiki.fysik.dtu.dk/ase/) — no need to write a `.py` file, just paste directly into your terminal.
- Commands are grouped by purpose below. Copy-paste as needed; most require no setup beyond having `ase` installed (`pip install ase`) and, for SLURM/module commands, being logged into your HPC cluster.

---

## 1. Format Conversion (ASE)

```bash
# Convert QE output (all relaxation steps) to a trajectory
python -O -c "from ase.io import read, write; write('output.xyz', read('input', format='espresso-out', index=':'))"

# Convert QE output (final structure only) to xyz/cif/vasp/etc.
python -O -c "from ase.io import read, write; write('output.vasp', read('input', format='espresso-out', index='-1'))"

# Convert a structure file (cif/xyz/vasp/POSCAR) into a QE input template
python -O -c "from ase.io import read, write; write('output.pwi', read('input'), format='espresso-in')"

# Generic format-agnostic conversion (ASE guesses from extension)
python -O -c "from ase.io import read, write; write('output', read('input'))"
```

## 2. Extracting Info from QE Runs (ASE)

```bash
# Print final total energy (eV)
python -O -c "from ase.io import read; a=read('input', format='espresso-out', index='-1'); print(a.get_potential_energy())"

# Print energy at every ionic step (convergence check)
python -O -c "from ase.io import read; [print(a.get_potential_energy()) for a in read('input', format='espresso-out', index=':')]"

# Print max force (eV/Å) on final structure
python -O -c "from ase.io import read; a=read('input', format='espresso-out', index='-1'); print(a.get_forces().max())"

# Print cell parameters of final structure
python -O -c "from ase.io import read; print(read('input', format='espresso-out', index='-1').cell.cellpar())"

# Count ionic steps in a relaxation
python -O -c "from ase.io import read; print(len(read('input', format='espresso-out', index=':')))"
```

## 3. Structure Inspection / Editing (ASE)

```bash
# Print chemical formula and number of atoms
python -O -c "from ase.io import read; a=read('input'); print(a.get_chemical_formula(), len(a))"

# Print atomic symbols + positions
python -O -c "from ase.io import read; a=read('input'); [print(s,p) for s,p in zip(a.get_chemical_symbols(), a.get_positions())]"

# Check magnetic moments (useful for nspin=2 setups)
python -O -c "from ase.io import read; a=read('input', format='espresso-out', index='-1'); print(a.get_magnetic_moments())"

# Center + wrap atoms into cell (fix PBC artifacts before viewing)
python -O -c "from ase.io import read, write; a=read('input'); a.wrap(); write('output', a)"

# Sort atoms by species (sometimes needed for clean QE input)
python -O -c "from ase.io import read, write; a=read('input'); write('output', a[a.numbers.argsort()])"
```

## 4. Visualization (ASE)

```bash
# Quick view with ASE GUI (needs a display)
python -O -c "from ase.io import read; from ase.visualize import view; view(read('input'))"

# Render a static image (png) instead of GUI, for headless HPC nodes
python -O -c "from ase.io import read, write; write('output.png', read('input'), rotation='10x,10y,10z')"
```

## 5. SLURM Job Monitoring

```bash
# Check your queued/running jobs
squeue -u $USER

# Check status of a specific job
squeue -j <jobid>

# Watch a job's queue status live (refresh every 5s)
watch -n 5 "squeue -u $USER"

# Cancel a job
scancel <jobid>

# Show detailed accounting info (runtime, memory, exit code) after a job finishes
sacct -j <jobid> --format=JobID,JobName,Elapsed,MaxRSS,State,ExitCode

# Tail a running job's output live
tail -f slurm-<jobid>.out
```

## 6. Quick QE Convergence/Status Checks (grep, no Python)

```bash
# Check if a QE run finished successfully
grep "JOB DONE" input

# Watch total energy converge (prints each SCF step's final energy line)
grep "!" input | tail -n 20

# Watch it live as the job progresses
tail -f input | grep --line-buffered "!"

# Grep max force per step (quick convergence check)
grep "Total force" input

# Check how many SCF iterations each step took (spot slow convergence)
grep -c "iteration #" input

# Pull wall time used so far
grep "WALL" input | tail -n 1

# Check for common QE crash signatures
grep -iE "error|not converged|crash" input
```

## 7. File / Disk Management on Shared HPC Storage

```bash
# Check disk usage of your project folder, sorted, human-readable
du -sh */ | sort -rh

# Check your quota (varies by cluster, common examples)
quota -s
lfs quota -u $USER /path/to/scratch

# Find all QE output files modified in the last day (e.g. just-finished jobs)
find . -name "*.out" -mtime -1

# Bulk copy relaxed structures out of many subfolders
find . -name "relax.out" -exec cp {} ./collected/{}_copy \;

# Compress a batch of finished job folders before archiving
tar -czf dopant_batch1.tar.gz job_*/

# Sync results back to local machine (run from local terminal, not HPC)
rsync -avz user@cluster:/path/to/results/ ./local_results/
```

## 8. Batch-Processing Many Structures at Once

```bash
# Loop over all QE outputs in subfolders and dump energies to a CSV
for d in job_*/; do e=$(python -O -c "from ase.io import read; print(read('$d/relax.out', format='espresso-out', index='-1').get_potential_energy())"); echo "$d,$e"; done > energies.csv

# Same, but flag any job that hasn't converged yet
for d in job_*/; do grep -q "JOB DONE" "$d/relax.out" || echo "$d NOT DONE"; done
```

## 9. Environment & Modules

```bash
# List loaded modules
module list

# Search for QE-related modules
module avail | grep -i espresso

# Load QE and see which binaries are exposed
module load quantum-espresso && which pw.x

# Check MPI ranks available before submitting
echo $SLURM_NTASKS
```

## 10. Resource Monitoring on the Node

```bash
# Live memory/CPU usage
htop        # or: top

# Free memory in human-readable form
free -h

# Check GPU usage (if your cluster uses GPU-accelerated QE)
nvidia-smi

# Keep a job running after you disconnect
tmux new -s dft_session
# reattach later with: tmux attach -t dft_session
```

## 11. Pseudopotential / Input Sanity Checks

```bash
# Confirm a pseudopotential file exists before submitting
ls $PSEUDO_DIR | grep -i <element>

# Quickly diff two QE input files (e.g. checking k-point changes between runs)
diff job1/input.pwi job2/input.pwi
```

---

### Notes

- `index=':'` gives every ionic step (useful for building training trajectories, e.g. for an active-learning ML potential); `index='-1'` gives only the final relaxed structure.
- `format='espresso-out'` is only needed when the filename doesn't end in something ASE auto-recognizes (e.g. `.pwo`) — if your log files are named like `relax.out`, specify the format explicitly.
- The batch-processing loops in Section 8 assume a folder structure like `job_1/relax.out`, `job_2/relax.out`, etc. — adjust the glob pattern (`job_*/`) to match your actual naming.
