"""Generate SLURM scripts for the bare (non-adsorbed) cluster runs.

Writes a ``run.sh`` into every ``<Metal>/{Tetra,Rhom}`` directory that contains
an ``INPUT_*.in``, defaulting to the newest input present.
"""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
METALS_DIR = PROJECT_ROOT / "05_metals"


def generate_run_script(tetra_dir):
    input_files = sorted(
        path for path in tetra_dir.glob("INPUT_*.in")
        if "adsorbed" not in path.name.lower()
    )
    if not input_files:
        print(f"Skipping {tetra_dir}: no cluster INPUT_*.in file found")
        return

    default_input = input_files[-1].name
    output_name = default_input.replace("INPUT_", "OUTPUT_").replace(".in", ".out")
    job_name = f"{tetra_dir.parent.name}_Tetra"

    content = """#!/bin/bash
#SBATCH --job-name='{job_name}'
#SBATCH --partition=high
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --requeue

set -euo pipefail

INPUT_FILE="${{1:-{default_input}}}"
if [[ -n "${{2:-}}" ]]; then
    OUTPUT_FILE="$2"
else
    OUTPUT_FILE="${{INPUT_FILE/INPUT_/OUTPUT_}}"
    OUTPUT_FILE="${{OUTPUT_FILE%.in}}.out"
fi

module load qe/7.2
mpiexec -np 8 pw.x -i "$INPUT_FILE" > "$OUTPUT_FILE"
""".format(job_name=job_name, default_input=default_input)

    run_path = tetra_dir / "run.sh"
    run_path.write_text(content)
    run_path.chmod(0o755)
    print(f"Generated {run_path.relative_to(PROJECT_ROOT)}")


def main():
    if not METALS_DIR.exists():
        raise FileNotFoundError(f"Metals directory missing: {METALS_DIR}")

    cluster_dirs = sorted(
        path
        for shape_name in ("Tetra", "Rhom")
        for path in METALS_DIR.glob(f"*/{shape_name}")
        if path.is_dir()
    )
    for cluster_dir in cluster_dirs:
        generate_run_script(cluster_dir)


if __name__ == "__main__":
    main()
