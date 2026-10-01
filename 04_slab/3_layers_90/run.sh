#!/bin/bash
#SBATCH --job-name='3x2_3_layers_90'
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

# FIXED: Removed '-i' and replaced it with '<' for proper input redirection
mpirun --mca orte_base_help_aggregate 0 -np 16 pw.x < 3x2_3_layers_90.in > 3x2_3_layers_90.out
