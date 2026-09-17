#!/bin/bash
#SBATCH --job-name='Pt4_Tetra'
#SBATCH --partition=high
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --requeue

set -euo pipefail

INPUT_FILE="${1:-INPUT_Pt4_rhom_ecut40.4_500.in}"
if [[ -n "${2:-}" ]]; then
    OUTPUT_FILE="$2"
else
    OUTPUT_FILE="${INPUT_FILE/INPUT_/OUTPUT_}"
    OUTPUT_FILE="${OUTPUT_FILE%.in}.out"
fi

module load qe/7.2
mpiexec -np 8 pw.x -i "$INPUT_FILE" > "$OUTPUT_FILE"
