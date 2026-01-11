#!/bin/bash
set -euo pipefail  # Go out if there is an error 

# The full code is already compiled in build,
# this bash only recompile if new modifications were made:
echo "COMPILING build: CPU in Leonardo:"
# Load modules 
module purge
module load profile/base
module load openmpi/4.1.6--gcc--12.2.0
module load cmake/

# Path to the build folder
BUILD_DIR="/leonardo_work/Sis25_baroni_0/gparedes/jacobi_perf/build"

cd "$BUILD_DIR"
echo "Compiling on host: $(hostname)"

# This part is not reaaly neeeded, because make check if something change or not:
#echo "[INFO-1/3] Compilin include ..." 
#make -C include/

echo "[INFO-2/3] Compiling jacobi CPU ..."
#make -j$(nproc) # Uses nproc processors
make -j$1

echo "[OK-3/3] Compilation  CPU leo done."

