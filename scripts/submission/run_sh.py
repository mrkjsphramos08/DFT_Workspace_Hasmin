"""Generate SLURM scripts for the adsorbed cluster runs.

Walks ``<Metal>/Tetra/Adsorbed/<profile>/<site>/`` and writes a ``run.sh`` for
each site that has an input but no output yet, so re-running is idempotent.
"""
import argparse
from pathlib import Path

from adsorption_config import DEFAULT_SLAB, SLAB_CHOICES, get_slab_profile
from adsorption_sites import VALID_SITES

# -------------------------------------------------------------------------
# Configuration & Path Resolution
# -------------------------------------------------------------------------
SCRIPT_PATH = Path(__file__).resolve()

PROJECT_ROOT = next(
    (p for p in [SCRIPT_PATH.parent] + list(SCRIPT_PATH.parents) if (p / "05_metals").exists()),
    SCRIPT_PATH.parent
)
METALS_DIR = PROJECT_ROOT / "05_metals"
 
def generate_slurm_script(input_in_path):
    """Generates run.sh in the same directory as the target INPUT_*.in file."""
    folder = input_in_path.parent
    site_name = folder.name  # e.g., 'centered', 'li', 'co', 'lico'
    
    # Locate 'Tetra' and Metal parent directory names dynamically
    tetra_dir = next((p for p in input_in_path.parents if p.name == "Tetra"), None)
    metal_name = tetra_dir.parent.name if tetra_dir else "QE"
    shape_name = tetra_dir.name if tetra_dir else "Job"
    
    # Format SLURM job name: e.g., Co4_Tetra_adsorbed_centered
    job_name = f"{metal_name}_{shape_name}_adsorbed_{site_name}"
    
    input_file_name = input_in_path.name
    output_file_name = input_file_name.replace("INPUT_", "OUTPUT_").replace(".in", ".out")
    
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name='{job_name}'
#SBATCH --partition=high
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --requeue

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

module load qe/7.2

mpiexec -np 16 pw.x -i {input_file_name} > {output_file_name}
"""
    
    # Derive a unique run.sh name from the template stem so multiple
    # input files (one per template) in the same site directory don't
    # overwrite each other's SLURM script.
    stem = input_in_path.name
    site = site_name
    marker = f"_{site}_"
    idx = stem.rfind(marker)
    if idx >= 0:
        template_stem = stem[idx + len(marker):].replace(".in", "")
    else:
        template_stem = stem.replace("INPUT_", "").replace(".in", "")
    run_sh_path = folder / f"run_{template_stem}.sh"
    with open(run_sh_path, "w") as f:
        f.write(slurm_content)

    print(f"Generated SLURM script: {run_sh_path.relative_to(PROJECT_ROOT)}")

def main():
    parser = argparse.ArgumentParser(description="Generate SLURM scripts for adsorbed-cluster QE inputs.")
    parser.add_argument(
        "--slab",
        choices=SLAB_CHOICES,
        default=DEFAULT_SLAB,
        help=f"Slab profile to process (default: {DEFAULT_SLAB})",
    )
    args = parser.parse_args()

    if not METALS_DIR.exists():
        print(f"Error: Directory '{METALS_DIR}' does not exist.")
        return

    output_name = get_slab_profile(args.slab)["output_name"]
    # Find all generated INPUT_*.in files in target subfolders.
    # Skip any site folder that already has a QE output file, since the job has already been run.
    input_files = []
    for input_path in METALS_DIR.rglob(f"Tetra/Adsorbed/{output_name}/*/*.in"):
        if input_path.parent.name not in VALID_SITES or not input_path.name.startswith("INPUT_"):
            continue

        # The slurm script redirects stdout to OUTPUT_<input_stem>.out;
        # check for that specific file so each variant is tracked independently.
        expected_output = input_path.name.replace("INPUT_", "OUTPUT_").replace(".in", ".out")
        has_output = (input_path.parent / expected_output).is_file()
        if not has_output:
            input_files.append(input_path)

    if not input_files:
        print("No fresh INPUT_*.in files found to create run.sh scripts for.")
        return

    for input_in_path in input_files:
        generate_slurm_script(input_in_path)

if __name__ == "__main__":
    main()
