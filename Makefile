# Compiler settings
CXX = g++
MPICXX = mpicxx
MPIRUN = mpirun

# execution file
TARGETS = MM_naive MM_strassen
OPENMP = naive_OpenMP strassen_OpenMP
MPI = naive_MPI strassen_MPI

# Default matrix size
R1 ?= 500
C1 ?= 500
R2 ?= 500
C2 ?= 500
THREADS ?= 4
PROCESS ?= 4

all: $(TARGETS)
openmp: $(OPENMP)
mpi: $(MPI)

MM_naive: MM_naive.cpp
	$(CXX) -o MM_naive MM_naive.cpp

MM_strassen: MM_strassen.cpp
	$(CXX) -o MM_strassen MM_strassen.cpp

naive_OpenMP: naive_OpenMP.cpp
	$(CXX) -o naive_OpenMP -fopenmp naive_OpenMP.cpp 

strassen_OpenMP: strassen_OpenMP.cpp
	$(CXX) -o strassen_OpenMP -fopenmp strassen_OpenMP.cpp 

naive_MPI: naive_MPI.cpp
	$(MPICXX) -o naive_MPI naive_MPI.cpp

strassen_MPI: strassen_MPI.cpp
	$(MPICXX) -o strassen_MPI strassen_MPI.cpp

run: all
	@echo "--- Running Naive ---"
	./MM_naive
	@echo "\n--- Running Strassen ---"
	./MM_strassen

runopenmp: openmp
	@echo "--- Running Naive ---"
	OMP_NUM_THREADS=$(THREADS) ./naive_OpenMP
	@echo "\n--- Running Strassen ---"
	OMP_NUM_THREADS=$(THREADS) ./strassen_OpenMP

runmpi: mpi
	@echo "--- Running Naive ---"
	$(MPIRUN) -np $(PROCESS) ./naive_MPI
	@echo "\n--- Running Strassen ---"
	$(MPIRUN) -np $(PROCESS) ./strassen_MPI

gen:
	python3 gen_matrix.py $(R1) $(C1) $(R2) $(C2)

clean:
	rm -f $(TARGETS) $(OPENMP) $(MPI) result_*.txt
