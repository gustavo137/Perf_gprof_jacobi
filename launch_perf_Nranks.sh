#!/bin/bash
set -euo pipefail

# We must run under srun
if [ -z "${SLURM_PROCID:-}" ]; then
  echo "SLURM_PROCID not set; run with srun."
  exit 1
fi

EXE="${EXE:-./build/bin/jacobi.x}"

# Frequency of sampling
# Default: 900 Hz
FREQ="${PERF_FREQ:-900}"

# Out directory per job (uses SCRATCH to avoid filling home/work)
OUTDIR="${PERF_OUTDIR:-$SCRATCH/jacobi/perf_${SLURM_JOB_ID}}"
mkdir -p "$OUTDIR"
umask 007

OUTFILE="${OUTDIR}/perf.data.rank${SLURM_PROCID}"

common_args=(
  -F "$FREQ"
  --output="$OUTFILE"
  --call-graph dwarf
  --delay 10
  -e cycles:u,instructions:u
  -g
)

# Forward all extra args (e.g., matrix size) to the executable
exec perf record "${common_args[@]}" -- "$EXE" "$@"
