#!/bin/bash
#SBATCH --job-name='Ni1_adsorbed_O'
#SBATCH --partition=high
#SBATCH --ntasks=40
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --exclude=n004

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

module load qe/7.2

mpiexec -np 40 pw.x -i Ni1_O.in > Ni1_O.out
