# Parallel Matrix Multiplication Assignment

Reference implementations of dense matrix multiplication used throughout the Parallel Computing assignment. The repository contains:

- Straightforward $O(n^3)$ multiplication for correctness and baseline timing.
- Strassen’s divide-and-conquer algorithm (serial, OpenMP, and MPI flavors) with automatic padding.
- OpenMP and MPI upgrades for the naive algorithm so you can study scalability on shared-memory and distributed clusters.
- Utility scripts to generate input matrices and log results for later analysis.

---

## � Repository layout

| Path | Description |
| --- | --- |
| `MM_naive.cpp` | Serial triple-loop implementation; writes `result_naive.txt`. |
| `MM_strassen.cpp` | Serial Strassen with power-of-two padding and configurable `LEAF_SIZE`; writes `result_strassen.txt`. |
| `naive_OpenMP.cpp`, `strassen_OpenMP.cpp` | Shared-memory versions using OpenMP; report thread count, execution time, and checksum to `result_openmp.txt`. |
| `naive_MPI.cpp`, `strassen_MPI.cpp` | Distributed-memory versions using MPI; `strassen_MPI` performs fully MPI-driven Strassen recursion and appends logs to `result_mpi.txt`. |
| `strassen_OpenMPI.cpp` | Spare sandbox for experimenting with hybrid ideas (not wired into the Makefile). |
| `gen_matrix.py` | Matrix generator (`python3 gen_matrix.py r1 c1 r2 c2`). |
| `matrix1`, `matrix2` | Sample inputs in the project’s plain-text format. |
| `Makefile` | Builds all variants and exposes convenience targets (`run`, `runopenmp`, `runmpi`, `gen`, `clean`). |
| `Assignment1_v1-1_031025.pdf` | Official assignment specification for reference. |
| `result_*.txt` | Timing/checksum logs produced by each variant (safe to delete). |

---

## ⚙️ Prerequisites

- Linux or another POSIX-like environment with `bash` and `make`.
- A C++17-capable compiler (`g++`, `clang++`, or similar) with OpenMP support.
- MPI toolchain (`mpicxx`, `mpirun`).
- Python 3 for matrix generation (optional but recommended).

Installations on Ubuntu (example):

```bash
sudo apt install build-essential libomp-dev openmpi-bin openmpi-common python3
```

---

## 🧱 Building

| Target | Command | Notes |
| --- | --- | --- |
| Serial baselines | `make` | Produces `MM_naive` and `MM_strassen`. |
| OpenMP variants | `make openmp` | Builds `naive_OpenMP` and `strassen_OpenMP`. |
| MPI variants | `make mpi` | Builds `naive_MPI` and `strassen_MPI`. |
| Everything | `make all openmp mpi` | Useful for cluster runs. |

Environment variables you can override per invocation:

- `THREADS` – Number of OpenMP threads used by `make runopenmp` (default `4`).
- `PROCESS` – Number of MPI ranks launched by `make runmpi` (default `4`).

---

## ▶️ Running the executables

### Serial reference

```bash
make run          # runs MM_naive then MM_strassen using matrix1/matrix2
```

Each executable prints matrix dimensions, execution time, and checksum, then writes the resulting matrix to `result_naive.txt` or `result_strassen.txt`.

### OpenMP

```bash
THREADS=8 make runopenmp
```

`OMP_NUM_THREADS` is injected automatically. Individual runs can also be launched manually:

```bash
export OMP_NUM_THREADS=8
./naive_OpenMP
./strassen_OpenMP
```

### MPI

```bash
PROCESS=6 make runmpi
```

Behind the scenes this calls:

```bash
mpirun -np 6 ./naive_MPI
mpirun -np 6 ./strassen_MPI
```

The MPI versions broadcast matrices from rank 0, divide work among ranks, gather the result, and append timing/checksum blocks to `result_mpi.txt`. The Strassen MPI version currently distributes the top recursion level, falling back to local Strassen for smaller blocks.

---

## 🗃️ Input format & generation

All programs expect matrices in plain text:

```
<rows> <cols>
v11 v12 ... v1c
...
vr1 vr2 ... vrc
```

Sample creation commands:

```bash
# create random 2048x2048 matrices (integers 0..9)
python3 gen_matrix.py 2048 2048 2048 2048

# Or via make (defaults to 500x500, override as needed)
R1=1024 C1=1024 R2=1024 C2=1024 make gen
```

Ensure the inner dimensions agree (`matrix1` columns == `matrix2` rows). Otherwise the programs exit with a descriptive error message.

---

## 🧮 Algorithms at a glance

| Variant | Highlights |
| --- | --- |
| Naive (serial/OpenMP/MPI) | Triple loop, no padding, deterministic output, easiest for correctness checks. |
| Strassen (serial) | Pads to the next power of two, recurses until `LEAF_SIZE`, then falls back to naive multiplication. |
| Strassen OpenMP | Creates OpenMP tasks for each of the seven Strassen sub-products and collapses them after `#pragma omp taskwait`. |
| Strassen MPI | Rank 0 prepares the seven sub-problems, distributes them across ranks, and merges the returned quadrants. A configurable `DISTRIBUTED_THRESHOLD` keeps tiny problems local. |

Checksums (sum of all matrix elements) are printed to help confirm that different variants yield identical outputs.

---

## 📊 Benchmarking & correctness checklist

1. **Warm-up**: Run the serial naive version to record a baseline time and checksum.
2. **Correctness**: Compare the checksum (or the entire `result_*.txt`) produced by each optimized variant against the serial baseline. Any mismatch indicates either data race issues or inconsistent padding.
3. **Scaling study**: Vary `THREADS` for OpenMP and `PROCESS` for MPI. Record (matrix size, processes/threads, runtime) triples in a spreadsheet or CSV.
4. **Matrix size sweep**: Use the generator to evaluate 100×100, 1 000×1 000, up to 10 000×10 000 matrices as requested in the assignment. Expect Strassen to pull ahead only after a certain crossover size.
5. **Library comparisons**: Optionally compare against BLAS (e.g., OpenBLAS or cuBLAS) to satisfy the “compare with existing library” requirement.

Suggested table template:

| Size | Variant | Threads/Processes | Time (s) | Speedup vs. serial |
| --- | --- | --- | --- | --- |
| 1024² | Naive OpenMP | 8 threads | 3.21 | 2.8× |

---

## 🧹 Maintenance commands

```bash
make clean          # remove all executables and result_*.txt
rm result_*.txt     # selectively clear logs
```

---

## 🛠️ Troubleshooting

- **File not found**: Verify `matrix1`/`matrix2` live in the project root and that your `mpirun` working directory is correct (use `mpirun -wd $(pwd)` if needed).
- **Dimension mismatch**: Check the first line of both inputs; the number of columns in the first matrix must equal the number of rows in the second.
- **OpenMP not parallelizing**: Make sure you compiled with `-fopenmp` (use the Makefile targets) and that `OMP_NUM_THREADS` is set to >1.
- **MPI ranks hang on exit**: Always allow rank 0 to reach the cleanup code. Killing a run abruptly may leave stray ranks; use `mpirun --timeout` or `pkill -f mpirun` as a last resort.
- **Strassen slow on tiny sizes**: Decrease `LEAF_SIZE`/`DISTRIBUTED_THRESHOLD`, or simply run the naive algorithm for small matrices.

---

## 🚀 Possible extensions

- Hybrid MPI + OpenMP (nested parallelism) to meet the “marriage” requirement from the assignment PDF.
- GPU kernels (CUDA, OpenCL, HIP) for extra credit.
- Automated CSV/JSON report generation and plotting scripts.
- Error-bound checks (e.g., relative Frobenius norm) when experimenting with floating-point reductions.

Happy benchmarking!
