# Matrix Multiplication Assignment

This assignment compares a baseline (\*O(n^3)\*) matrix multiplication with an optimized Strassen-based approach. Both programs are implemented in C++ and share the same input/output format so you can benchmark accuracy and runtime side-by-side.

## 📁 Project layout

| Path | Purpose |
| --- | --- |
| `MM_naive.cpp` | Classic triple-loop multiplication reading matrices from disk and writing `result_naive.txt` |
| `MM_strassen.cpp` | Strassen implementation with power-of-two padding and a leaf-size cutoff for small subproblems; writes `result_strassen.txt` |
| `matrix1`, `matrix2` | Sample matrices stored as plain text (first line: rows cols, followed by row-major elements) |
| `Makefile` | Builds both binaries and provides convenience targets for running/cleaning |

## ⚙️ Prerequisites

- A C++17-compatible compiler (tested with `g++` on Linux)
- `make`

## 🚧 Build

```bash
make
```

This compiles two executables in the project root: `MM_naive` and `MM_strassen`.

## ▶️ Run

Use the provided `run` target to execute both programs sequentially with the current `matrix1`/`matrix2` data:

```bash
make run
```

Or run either binary manually:

```bash
./MM_naive
./MM_strassen
```

Each program prints the detected matrix dimensions, starts a high-resolution timer, reports the elapsed time in seconds, and writes the output product to `result_naive.txt` or `result_strassen.txt` (first line contains the result dimensions).

## 📦 Input format

1. Plain-text file.
2. First line: `<rows> <cols>`.
3. Remaining lines: row-major elements separated by spaces/newlines.

Example (`matrix1`):

```
1 3
3 4 2
```

`matrix2` follows the same convention with its own dimensions and values. The inner dimensions must match (`c1 == r2`); otherwise the programs will terminate early with an error message.

## 🧮 Algorithms

### Naive multiplication
- Uses the classic triple-nested loop.
- Works for all rectangular matrix shapes without padding.

### Strassen multiplication
- Automatically pads both input matrices to the next power of two so recursion works on square blocks.
- Recursion stops at a configurable `LEAF_SIZE` (currently 64), then falls back to a naive multiplication for the base case.
- After computation, the padded matrix is trimmed back to the original output dimensions.

## 📊 Benchmarking tips

- Replace `matrix1` and `matrix2` with larger datasets following the same format to stress-test both implementations.
- Run multiple iterations and average the reported execution times for more stable measurements.
- Because `MM_strassen` introduces recursion and temporary allocations, the crossover point where it outperforms the naive version will depend on matrix size and hardware cache behavior.

## 🧹 Housekeeping

```bash
make clean
```

Removes both executables and any previously generated `result_*.txt` files.

## 🛠️ Troubleshooting

- **Cannot open file**: Ensure `matrix1` and `matrix2` exist in the project root and the process has read permissions.
- **Dimension mismatch**: Check the first line of each file; the number of columns in `matrix1` must equal the number of rows in `matrix2`.
- **Slow Strassen run on small matrices**: Lower the `LEAF_SIZE` constant or switch to the naive algorithm for small test cases.

## ✅ Next steps

- Add additional datasets (random, dense, sparse) to explore performance trends.
- Extend the programs to emit CSV timing summaries for automated benchmarking.
- Parallelize either algorithm (OpenMP, MPI, CUDA) for larger assignments.
