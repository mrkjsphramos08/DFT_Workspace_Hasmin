#!/bin/bash
#SBATCH --job-name='Pt4_Tetra_adsorbed_lili'
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

mpiexec -np 16 pw.x -i INPUT_adsorbed_Pt4_Tetra_latest_lili.in > OUTPUT_adsorbed_Pt4_Tetra_latest_lili.out
