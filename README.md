# Profiliong Instructions

## Normal compilation
```
module laod cmake/
module load openmpi/4.1.6--gcc--12.2.0

mkdir -p build
cd build
cmake ..
make -j
```

## Using gprof 

```
cd build
cmake -DENABLE_GPROF=ON ..
make -j
```
Clean the previous gmon.out files if any:
```
rm -rf gmon.out gmon.out.*
```
### Runing with one process MPI (1 Rank)
We run the program with:
```
srun -n 1 ./build/bin/jacobi.x 1000 > out/jacobi_1N.dat
```
Then we can generate the gprof basic report with:
```
gprof ./build/bin/jacobi.x gmon.out > out/gprof_report.txt
```
### Generate the call graph with Graphviz(dot) / gprof2dot
We need to install graphviz and gprof2dot if not already installed:
```
pip install gprof2dot
```
Then we can generate the call graph with:
```
gprof -b ./build/bin/jacobi.x gmon.out | gprof2dot -f prof | dot -Tpdf -o callgraph.pdf
```
or for PNG format:
```
gprof -b ./build/bin/jacobi.x gmon.out | gprof2dot -f prof | dot -Tpng -o callgraph.png
```

### Running with N processes MPI



### Install graphviz/gprof2dot on cluster
After doing `modmap -m dot graphviz` we check that in leonardo is already installed:
```
which dot
```
To install gprof2dot we can use a virtual environment:
```
module load python/3.11.7
python -m venv ~/myenv
source ~/myenv/bin/activate
pip install --upgrade pip
pip install gprof2dot
python3 -m pip install --user gprof2dot
```
Then we need to add the local bin to our PATH:
```
export PATH="$HOME/.local/bin:$PATH"
```
Now you can use gprof2dot in this env. When done, you can exit the env with:
```
deactivate
```
Then we can use gprof2dot as shown above.
f you want to remove the env do: 
```
rm -rf ~/myenv
```
## Using perf 

```
module laod cmake/
module load openmpi/4.1.6--gcc--12.2.0

rm -rf build
cmake -S . -B build -DENABLE_PERF=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Generating flame graphs from perf data
We need to clone `git clone git@github.com:brendangregg/FlameGraph.git` to generate flame graphs from perf data and `export PATH="$PATH:$PWD/FlameGraph"`.
Then we can run the program with:
```
# Option 2: Using perf to collect performance data
# ----- Run with perf event one process MPI -----
# perf record --delay 10 -F 999 -g --call-graph dwarf -e cycles:u,instructions:u -- \
#   ./build/bin/jacobi.x 1000 > out/jacobi_1N.dat
```
to collect the perf data, (check the `job.slurm`).
To optain the graphs:
```
perf script -i perf.data | stackcollapse-perf.pl > out.folded
flamegraph.pl --title "Jacobi Serial (1 rank)" out.folded > flamegraph_serial.svg
```
directly in the terminal or run the `do_perf_graphs_serial.sh` to generate the flame graphs.
```
#!/bin/bash
set -euo pipefail  # Go out if there is an error
# Generate flamegraph for Jacobi Serial (1 rank)

echo "Generating flamegraph for Jacobi Serial (1 rank)..."

perf script -i perf.data | stackcollapse-perf.pl > out.folded
flamegraph.pl --title "Jacobi Serial (1 rank)" out.folded > flamegraph_serial.svg

echo "Flamegraph generated: flamegraph_serial.svg"

# Clean up intermediate files
rm out.folded perf.data
```


To open the `perf_flamegraph.svg` in a browser to see the flame graph we can use:
```
google-chrome flamegraph_serial.svg`           #on Linux
open -a "Brave Browser" flamegraph_serial.svg` # on macOS
```
we need to copy the `flamegraph_serial.svg` to our local machine if we are working on a remote server.

### Using perf with N MPI processes
Here is more convenient a wrapper script to run perf with N MPI processes `launch_perf_Nranks.sh`:
```
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

exec perf record "${common_args[@]}" -- "$EXE"
```
we give execute permission to the script `chmod +x launch_perf_Nranks.sh` and then we can use it in the `job.slurm` as follows:
```
# Option 3: Using perf with N processes MPI
# ----- Run with perf event N process MPI -----
srun -n 4 ./launch_perf_Nranks.sh 1000 > out/jacobi_1N.dat
```

Then the results are stored in the `SCRATCH/jacobi/perf_${SLURM_JOB_ID}` directory.
```
ls $SCRATCH/jacobi/perf_${SLURM_JOB_ID}/perf.data.rank*
```
#### Postprocessing perf data from N MPI processes
Option A: merge flamegraph from rank 0 ... 3
We can use parallel to decode several perf.data.rank* in parallel and merges them into a single .perf.
On an interactive node (recommended if heavy): 
```
srun -A ICT25_MHPC -p dcgp_usr_prod -N1 -n4 --cpus-per-task=1 --mem=0g --gres=tmpfs:300g --time=01:00:00 --pty bash
```
Then run:
```
OUTDIR=$SCRATCH/jacobi/perf_${SLURM_JOB_ID}
cd "$OUTDIR"

# Decode and merge all ranks (0 ,..., 3) in parallel
ls perf.data.rank{0..3} | parallel -j4 'perf script -i {} --max-stack 128 --no-demangle' > subset_0_3.perf

stackcollapse-perf.pl subset_0_3.perf > rank0_3.folded
flamegraph.pl --title "Jacobi MPI 4 tasks" rank0_3.folded > flamegraph_mpi_4tasks.svg
```
or directly run the job script `job-merge-perf-ranks.slurm`. 
