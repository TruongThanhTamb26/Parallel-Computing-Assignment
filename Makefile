# Compiler settings
CXX = g++
MPICXX = mpicxx
MPIRUN = mpirun

# execution file
NORMAL = MM_naive MM_strassen
OPENMP = naive_OpenMP strassen_OpenMP
MPI = naive_MPI strassen_MPI
OPENMPI = naive_OpenMPI strassen_OpenMPI

# Default matrix size
R1 ?= 500
C1 ?= 500
R2 ?= 500
C2 ?= 500
THREADS ?= 4
PROCESS ?= 4

normal: $(NORMAL)
openmp: $(OPENMP)
mpi: $(MPI)
openmpi: $(OPENMPI)
all: $(NORMAL) $(OPENMP) $(MPI) $(OPENMPI)

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

naive_OpenMPI: naive_OpenMPI.cpp
	$(MPICXX) -fopenmp -o naive_OpenMPI naive_OpenMPI.cpp

strassen_OpenMPI: strassen_OpenMPI.cpp
	$(MPICXX) -fopenmp -o strassen_OpenMPI strassen_OpenMPI.cpp

runnormal: normal
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

runopenmpi: openmpi
	@echo "--- Running Naive ---"
	OMP_NUM_THREADS=$(THREADS) $(MPIRUN) -np $(PROCESS) ./naive_OpenMPI 
	@echo "\n--- Running Strassen ---"
	OMP_NUM_THREADS=$(THREADS) $(MPIRUN) -np $(PROCESS) ./strassen_OpenMPI

runall: runnormal runopenmp runmpi runopenmpi

#run with cluster
runmpicluster: mpi
	@echo "--- Running Naive ---"
	$(MPIRUN) -f hosts.txt -n $(PROCESS) ./naive_MPI
	@echo "\n--- Running Strassen ---"
	$(MPIRUN) -f hosts.txt -n $(PROCESS) ./strassen_MPI

runopenmpicluster: openmpi
	@echo "--- Running Naive ---"
	OMP_NUM_THREADS=$(THREADS) $(MPIRUN) -f hosts.txt -np $(PROCESS) ./naive_OpenMPI
	@echo "\n--- Running Strassen ---"
	OMP_NUM_THREADS=$(THREADS) $(MPIRUN) -f hosts.txt -np $(PROCESS) ./strassen_OpenMPI


gen:
	python3 gen_matrix.py $(R1) $(C1) $(R2) $(C2)

clean:
	rm -f $(NORMAL) $(OPENMP) $(MPI) $(OPENMPI) result_*.txt
