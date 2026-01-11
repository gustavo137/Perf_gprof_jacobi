#!/bin/bash
set -euo pipefail  # Go out if there is an error
# Generate flamegraph for Jacobi Serial (1 rank)

echo "Generating flamegraph for Jacobi Serial (1 rank)..."
export PATH="$PATH:$PWD/FlameGraph"
# export PATH="$PATH:$SLURM_SUBMIT_DIR/FlameGraph"

perf script -i perf.data | stackcollapse-perf.pl > out.folded
flamegraph.pl --title "Jacobi Serial (1 rank)" out.folded > flamegraph_serial.svg

echo "Flamegraph generated: flamegraph_serial.svg"

# Clean up intermediate files
rm out.folded perf.data
