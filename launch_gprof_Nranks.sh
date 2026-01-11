#!/bin/bash
set -euo pipefail

: "${SLURM_PROCID:?Run under srun}"

EXE="${EXE:-./build/bin/jacobi.x}"

OUTDIR="${GPROF_OUTDIR:-$SCRATCH/jacobi/gprof_${SLURM_JOB_ID}}"
mkdir -p "$OUTDIR"
umask 007

# gprof will create files like: ${GMON_OUT_PREFIX}.<pid>
export GMON_OUT_PREFIX="${OUTDIR}/gmon.rank${SLURM_PROCID}"

exec "$EXE" "$@"

