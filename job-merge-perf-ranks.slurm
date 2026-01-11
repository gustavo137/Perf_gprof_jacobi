#!/bin/bash
#SBATCH --job-name=jacobi1N              # Job name
##SBATCH --output=jacobi_1N.datg
#SBATCH --error=out/job.err            # std-error file
#SBATCH --output=out/job.out           # std-output file
#SBATCH -N 1                             # Number of nodes
#SBATCH --ntasks-per-node=32             # CPU MPI
#SBATCH --cpus-per-task=1                # CPU OpenMP
#SBATCH --time=00:30:00                  # Time limit hrs:min:sec
#SBATCH -A CMPNS_sissabar
#SBATCH -p boost_usr_prod

# echo "Running jacobi"

# Load the required modules
module purge
module load openmpi/4.1.6--gcc--12.2.0
#module load cmake/

export DIR_NAME=${SLURM_JOB_NAME}
mkdir -p out || true

# export PATH="$PATH:$PWD/FlameGraph"
export PATH="$PATH:$SLURM_SUBMIT_DIR/FlameGraph"

# Export the path to the executable
export EXE=$SLURM_SUBMIT_DIR/build/bin/jacobi.x

## Control of threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Perf settings
export PERF_OUTDIR=$SCRATCH/jacobi/perf_${SLURM_JOB_ID}

##for timing slurm
start_time=$SECONDS

# Option 3: Using perf with N processes MPI
# ----- Run with perf event N process MPI -----
srun -n 4 ./launch_perf_Nranks.sh 1000 > out/jacobi_1N.dat
#
# To aboid stdout overhead, redirect output to a file
# srun -n 4 --output=out/jacobi_rank%t.out ./launch_perf_Nranks.sh 1000

# echo "Job Finished"
cd "$PERF_OUTDIR"

# Merge ranks 0..3 without GNU parallel
rm -f subset_0_3.perf
for r in 0 1 2 3; do
  perf script -i "perf.data.rank${r}" --max-stack 128 --no-demangle >> subset_0_3.perf
done

stackcollapse-perf.pl subset_0_3.perf > rank0_3.folded
flamegraph.pl --title "Jacobi MPI 4 tasks" rank0_3.folded > "$SLURM_SUBMIT_DIR/flamegraph_mpi_4tasks.svg"

rm -rf "$PERF_OUTDIR"

end_time=$SECONDS
runtime=$((end_time - start_time))
echo "Run time: ${runtime} seconds ($(awk "BEGIN {printf \"%.2f\", ${runtime}/60}") minutes)." \
  > "$SLURM_SUBMIT_DIR/out/time_run_${SLURM_JOB_ID}.txt"

cd "$SLURM_SUBMIT_DIR"
mv out "out_${DIR_NAME}"

