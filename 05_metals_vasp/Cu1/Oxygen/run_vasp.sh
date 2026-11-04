#!/usr/bin/env bash
# Submit from the calculation folder with: sbatch run_vasp.sh
#SBATCH --job-name=Cu1_oxygen
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --output=Output.%j.out
#SBATCH --error=Error.%j.err

set -eo pipefail

if [[ -z ${SLURM_JOB_ID:-} || -z ${SLURM_SUBMIT_DIR:-} ]]; then
    echo 'Run this script through sbatch from the calculation directory.' >&2
    exit 2
fi
if [[ ${SLURM_NTASKS:-0} != 40 || ${SLURM_CPUS_PER_TASK:-1} != 1 || ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 ]]; then
    echo 'Expected one node, 40 MPI ranks, and one CPU per rank.' >&2
    exit 2
fi

cd -- "$SLURM_SUBMIT_DIR"
for input in INCAR POSCAR KPOINTS POTCAR; do
    if [[ ! -f $input || ! -r $input || ! -s $input ]]; then
        printf 'Missing, empty, or unreadable input: %s\n' "$input" >&2
        exit 2
    fi
done

oneapi_setvars=${ONEAPI_SETVARS:-/opt/intel/oneapi/setvars.sh}
if [[ ! -r $oneapi_setvars ]]; then
    printf 'Cannot read oneAPI environment script: %s\n' "$oneapi_setvars" >&2
    exit 2
fi

set +e
source "$oneapi_setvars"
setup_status=$?
set -e
if (( setup_status != 0 )); then
    printf 'oneAPI setup failed with status %s\n' "$setup_status" >&2
    exit "$setup_status"
fi
set -u

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export I_MPI_HYDRA_BOOTSTRAP=slurm

vasp_exe=${VASP_EXE:-/home/Software/VASP/vasp.6.5.1/bin/vasp_std}
command -v mpiexec >/dev/null || { echo 'Intel MPI mpiexec not found.' >&2; exit 2; }
command -v "$vasp_exe" >/dev/null || { echo 'VASP executable not found.' >&2; exit 2; }

printf 'Job: %s\nNodes: %s\nMPI ranks: %s\n' "$SLURM_JOB_ID" "${SLURM_JOB_NODELIST:-unknown}" "$SLURM_NTASKS"
printf 'MPI launcher: %s\nVASP executable: %s\n' "$(command -v mpiexec)" "$(command -v "$vasp_exe")"
printf 'Start: %s\n' "$(date -Is)"

status=0
mpiexec -np "$SLURM_NTASKS" "$vasp_exe" || status=$?
printf 'VASP exit code: %s\nFinished: %s\n' "$status" "$(date -Is)"
exit "$status"
