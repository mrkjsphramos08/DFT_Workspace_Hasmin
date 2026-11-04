#!/usr/bin/env bash
# Submit from the calculation folder with: sbatch run_vasp.sh
# Site assumption from the referenced chat: Hasmin's high partition.
#SBATCH --job-name=LCO_bulk
#SBATCH --partition=high
# One node must provide at least 40 allocated CPUs; no multi-node mapping here.
#SBATCH --nodes=1
# Forty MPI ranks implements the user's allocation; KPAR=4 gives 10 per group.
#SBATCH --ntasks=40
# One CPU per rank: use pure MPI, with one OpenMP/MKL thread per rank below.
#SBATCH --cpus-per-task=1
# Scheduler output filenames use %j (job ID); no shell variables in directives.
#SBATCH --output=Output.%j.out
#SBATCH --error=Error.%j.err
# Time, memory, account and QOS use site defaults; add required site directives
# ABOVE the first executable line. The actual cluster policy was not queried.

set -eo pipefail
# Refuse accidental execution on a login node or an allocation of another size.
if [[ -z ${SLURM_JOB_ID:-} || -z ${SLURM_SUBMIT_DIR:-} ]]; then
    echo 'Run this script through sbatch from the calculation directory.' >&2
    exit 2
fi
if [[ ${SLURM_NTASKS:-0} != 40 || ${SLURM_CPUS_PER_TASK:-1} != 1 || ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 ]]; then
    echo 'Expected one node, 40 MPI ranks, and one CPU per rank.' >&2
    exit 2
fi
# Work in the submit directory, where VASP reads its fixed input filenames.
cd -- "$SLURM_SUBMIT_DIR"
for input in INCAR POSCAR KPOINTS POTCAR; do
    if [[ ! -f $input || ! -r $input || ! -s $input ]]; then
        printf 'Missing, empty, or unreadable input: %s\n' "$input" >&2
        exit 2
    fi
done

# Intel oneAPI supplies runtime libraries and Intel MPI, not a VASP license
# or necessarily the VASP executable. Override only if the site uses another
# trusted oneAPI installation: export ONEAPI_SETVARS=/site/path/setvars.sh
oneapi_setvars=${ONEAPI_SETVARS:-/opt/intel/oneapi/setvars.sh}
if [[ ! -r $oneapi_setvars ]]; then
    printf 'Cannot read oneAPI environment script: %s\n' "$oneapi_setvars" >&2
    exit 2
fi
# Vendor setup scripts can reference unset variables; enable nounset afterward.
set +e
source "$oneapi_setvars"
setup_status=$?
set -e
if (( setup_status != 0 )); then
    printf 'oneAPI setup failed with status %s\n' "$setup_status" >&2
    exit "$setup_status"
fi
set -u

# Avoid 40 ranks each spawning additional CPU threads.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
# Ask Intel MPI Hydra to launch through Slurm, keeping scheduler allocation.
# Do not wrap this mpiexec command inside another srun command.
export I_MPI_HYDRA_BOOTSTRAP=slurm

# The standard executable is required for a Gamma-centered mesh with many
# k points. vasp_gam is a Gamma-only executable and is inappropriate here.
# A site may require loading its VASP module first; confirm ABI compatibility.
# Absolute path to the installed VASP 6.5.1 executable.
# Required because vasp_std is not currently in PATH.
vasp_exe=${VASP_EXE:-/home/Software/VASP/vasp.6.5.1/bin/vasp_std}
command -v mpiexec >/dev/null || { echo 'Intel MPI mpiexec not found.' >&2; exit 2; }
command -v "$vasp_exe" >/dev/null || { echo 'VASP executable not found.' >&2; exit 2; }
printf 'Job: %s\nNodes: %s\nMPI ranks: %s\n' "$SLURM_JOB_ID" "${SLURM_JOB_NODELIST:-unknown}" "$SLURM_NTASKS"
printf 'MPI launcher: %s\nVASP executable: %s\n' "$(command -v mpiexec)" "$(command -v "$vasp_exe")"
printf 'Start: %s\n' "$(date -Is)"

# Preserve the MPI/VASP exit code in the Slurm job status, even on failure.
# A zero exit code still requires inspecting electronic/ionic convergence.
status=0
mpiexec -np "$SLURM_NTASKS" "$vasp_exe" || status=$?
printf 'VASP exit code: %s\nFinished: %s\n' "$status" "$(date -Is)"
exit "$status"
