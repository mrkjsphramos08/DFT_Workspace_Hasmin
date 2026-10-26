#!/bin/bash
#SBATCH --job-name='Pt1_coco_441'
#SBATCH --partition=high
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --exclude=n004

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

module load qe/7.2

mpiexec -np 8 pw.x -i ../INPUT_adsorbed_Pt1_latest_coco_44.1_template.in > OUTPUT_adsorbed_Pt1_latest_coco_44.1_template.out
