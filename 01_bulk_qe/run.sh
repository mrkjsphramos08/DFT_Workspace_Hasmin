#!/bin/bash
#SBATCH --job-name=LCO_bulk
#SBATCH --partition=high
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mail-type=ALL

# Echo job details to standard output
echo "SLURM_JOBID = $SLURM_JOBID"
echo "SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "SLURM_NNODES = $SLURM_NNODES"
echo "Working directory = $SLURM_SUBMIT_DIR"

# Reset environment and load dependencies
module purge
module load openmpi
module load qe/7.2

# Prevent OpenMPI transport failure if InfiniBand is missing
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp

cd "$SLURM_SUBMIT_DIR"

# Execute pw.x via mpirun using the allocated task count
mpirun -np $SLURM_NTASKS pw.x -i LCO_bulk.in > LCO_bulk.out