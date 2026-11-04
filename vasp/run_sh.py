"""Generate SLURM scripts for VASP runs.

Mirrors ``run_sh.py`` / ``run_sh_pt.py`` (QE) but emits ``run_vasp``.sh for each
VASP input directory.  A run script is only written when the site directory
already contains an ``INCAR`` but no ``OSZIC``/``OUTCAR`` yet, so re-running is
idempotent and never clobbers an already-finished job.

POTCARs are deliberately **not** committed (licensed PBE PAW datasets, same
policy as the QE ``pseudo_dir`` on the cluster).  The generated script expects
``$VASP_PP_PATH`` to point at a PBE PAW library and assembles ``POTCAR`` from
the ``POTCAR_SPECIES`` file written by :mod:`vasp.generate_inputs`, before
launching ``vasp_std``.
"""
import argparse
from pathlib import Path

from lco.paths import VASP_METALS_DIR

# Task counts mirror the QE campaign (see README "SLURM resources" table).
NTASKS = {"cluster": 8, "adsorbed": 16, "slab": 32, "bulk": 32}

SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name='{job_name}'
#SBATCH --partition={partition}
#SBATCH --ntasks={ntasks}
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --requeue

set -euo pipefail

echo "SLURM_JOBID=$SLURM_JOBID"
echo "SLURM_JOB_NODELIST=$SLURM_JOB_NODELIST"
echo "SLURM_NNODES=$SLURM_NNODES"
echo "Working directory=$SLURM_SUBMIT_DIR"

module purge
module load openmpi
module load vasp/{version}          # <-- set your cluster's VASP module

# Assemble the (git-ignored) POTCAR from your PBE PAW library.
: "${{VASP_PP_PATH:?set VASP_PP_PATH to your PBE PAW library root}}"
if [[ ! -f POTCAR ]]; then
  mapfile -t SPECIES < <(awk '{{print $1}}' POTCAR_SPECIES)
  : > POTCAR
  for s in "${{SPECIES[@]}}"; do
    cat "$VASP_PP_PATH/$s/POTCAR" >> POTCAR   # adjust subdirectory (e.g. PBE/$s)
  done
fi

mpiexec -n $SLURM_NTASKS vasp_std > vasp.out
"""


def generate_slurm_script(vasp_dir, job_name=None, partition="high",
                          ntasks=None, version="6.4.2", overwrite=False):
    """Write ``run_vasp.sh`` into ``vasp_dir``; idempotent (skips if OSZIC/OUTCAR exist)."""
    vasp_dir = Path(vasp_dir)
    already_done = (vasp_dir / "OSZIC").exists() or (vasp_dir / "OUTCAR").exists() or \
        any(vasp_dir.glob("OSZIC*")) or any(vasp_dir.glob("OUTCAR*"))
    run_path = vasp_dir / "run_vasp.sh"
    if run_path.exists() and not overwrite and not already_done:
        return run_path

    if job_name is None:
        job_name = vasp_dir.parent.name
    if ntasks is None:
        ntasks = (NTASKS["adsorbed"] if "Adsorbed" in vasp_dir.parts
                  else NTASKS["cluster"] if any(p in vasp_dir.parts for p in ("Co4", "Pt4", "Cu4", "Ni4", "Pt1", "Ni1"))
                  else NTASKS["bulk"])

    content = SLURM_TEMPLATE.format(
        job_name=job_name, partition=partition, ntasks=ntasks, version=version
    )
    vasp_dir.mkdir(parents=True, exist_ok=True)
    run_path.write_text(content, encoding="utf-8")
    run_path.chmod(0o755)
    return run_path


def iter_vasp_dirs(vasp_root):
    """Yield directories that contain an INCAR under ``vasp_root``."""
    for incar in sorted(vasp_root.rglob("INCAR")):
        yield incar.parent


def main(argv=None):
    ap = argparse.ArgumentParser(description="Generate SLURM scripts for VASP runs.")
    ap.add_argument("--vasp-dir", default=str(VASP_METALS_DIR))
    ap.add_argument("--partition", default="high")
    ap.add_argument("--ntasks", type=int, default=None)
    ap.add_argument("--version", default="6.4.2", help="module name / version of VASP")
    ap.add_argument("--job-name", default=None)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args(argv)

    vasp_root = Path(args.vasp_dir)
    written = 0
    for d in iter_vasp_dirs(vasp_root):
        generate_slurm_script(d, job_name=args.job_name, partition=args.partition,
                              ntasks=args.ntasks, version=args.version,
                              overwrite=args.overwrite)
        written += 1
        print(f"  wrote run_vasp.sh in {d}")
    print(f"generated {written} run_vasp.sh scripts")


if __name__ == "__main__":
    main()
