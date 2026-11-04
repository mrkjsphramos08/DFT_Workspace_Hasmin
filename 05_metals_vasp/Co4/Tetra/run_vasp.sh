#!/bin/bash
#SBATCH --job-name='Co4_Tetra_vasp'
#SBATCH --partition=high
#SBATCH --ntasks=8
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
module load vasp/6.4.2          # <-- set your cluster's VASP module

# Assemble the (git-ignored) POTCAR from your PBE PAW library.
: "${VASP_PP_PATH:?set VASP_PP_PATH to your PBE PAW library root}"
if [[ ! -f POTCAR ]]; then
  mapfile -t SPECIES < <(awk '{print $1}' POTCAR_SPECIES)
  : > POTCAR
  for s in "${SPECIES[@]}"; do
    cat "$VASP_PP_PATH/$s/POTCAR" >> POTCAR   # adjust subdirectory (e.g. PBE/$s)
  done
fi

mpiexec -n $SLURM_NTASKS vasp_std > vasp.out
