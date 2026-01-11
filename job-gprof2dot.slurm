#!/bin/bash
#SBATCH --job-name=jacobi1N              # Job name
##SBATCH --output=jacobi_1N.datg
#SBATCH --error=out/job.err            # std-error file
#SBATCH --output=out/job.out           # std-output file
#SBATCH -N 1                             # Number of nodes
#SBATCH --ntasks-per-node=32             # CPU MPI
#SBATCH --cpus-per-task=1                # CPU OpenMP
#SBATCH --time=00:05:00                  # Time limit hrs:min:sec
#SBATCH -A CMPNS_sissabar
#SBATCH -p boost_usr_prod

# echo "Running jacobi"

# Load the required modules
module purge
module load openmpi/4.1.6--gcc--12.2.0
module load python/3.11.7
#module load cmake/

export DIR_NAME=${SLURM_JOB_NAME}
mkdir -p out || true

# Export the path to the executable
export EXE=$SLURM_SUBMIT_DIR/build/bin/jacobi.x

## Control of threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

#  ===============================================================
# -----  Run MPI program -----
# mpirun -np 32 "$exe" 1000 > out/jacobi_1N.dat
# srun --cpu_bind=cores "$EXE" 1000 > out/jacobi_1N.dat
#  ===============================================================

# Option 1: Using gprof to collect performance data
# ----- Run with gprof event one process MPI (serial) -----
# srun -n 1 "$EXE" 1000 > out/jacobi_1N.dat
# This will generate the gmon.out file in the working directory
# Now you can run `do_perf_graph_serial.sh` to generate the report after the job is done.
# ./do_perf_graph_serial.sh
#  ===============================================================
#
# Option 2: Using gprof to visualize data
OUTDIR=$SCRATCH/jacobi/gprof_${SLURM_JOB_ID}
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# ----- Run with gprof with N processes MPI -----
srun -n 4 ./launch_gprof_Nranks.sh 1000 > out/jacobi_1N.dat

cd "$OUTDIR"


for f in gmon.rank*; do
  gprof -b "$EXE" "$f" > "gprof.${f}.txt"
done

# Now generate call graphs in PDF format
source ~/env_gprof2dot/bin/activate

for f in gmon.rank*; do
  gprof -b "$EXE" "$f" \
    | gprof2dot -f prof \
    | dot -Tpdf -o "callgraphN.${f}.pdf"
done
deactivate
# Copy results back to submit dir (so you can download/inspect easily)
cp -a "$OUTDIR" "$SLURM_SUBMIT_DIR/gprof_${SLURM_JOB_ID}"
# Optional cleanup of scratch dir:
rm -rf "$OUTDIR"
cd "$SLURM_SUBMIT_DIR/gprof_${SLURM_JOB_ID}"
rm gmon.rank*
rm gprof.*.txt
#  ===============================================================
cd "$SLURM_SUBMIT_DIR"
# mv out  to out_"${DIR_NAME}"
mv out "out_${DIR_NAME}"

