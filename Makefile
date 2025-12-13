# Compiler settings
CXX = g++

# execution file
TARGETS = MM_naive MM_strassen
OPENMP = naive_OpenMP strassen_OpenMP

# Default matrix size
R1 ?= 500
C1 ?= 500
R2 ?= 500
C2 ?= 500
THREADS ?= 4

all: $(TARGETS)
openmp: $(OPENMP)

MM_naive: MM_naive.cpp
	$(CXX) -o MM_naive MM_naive.cpp

MM_strassen: MM_strassen.cpp
	$(CXX) -o MM_strassen MM_strassen.cpp

naive_OpenMP: OpenMP/naive.cpp
	$(CXX) -o OpenMP/naive -fopenmp OpenMP/naive.cpp 

strassen_OpenMP: OpenMP/strassen.cpp
	$(CXX) -o OpenMP/strassen -fopenmp OpenMP/strassen.cpp 

run: all
	@echo "--- Running Naive ---"
	./MM_naive
	@echo "\n--- Running Strassen ---"
	./MM_strassen

runopenmp: openmp
	@echo "--- Running Naive ---"
	OMP_NUM_THREADS=$(THREADS) ./OpenMP/naive
	@echo "\n--- Running Strassen ---"
	OMP_NUM_THREADS=$(THREADS) ./OpenMP/strassen

gen:
	python3 gen_matrix.py $(R1) $(C1) $(R2) $(C2)

clean:
	rm -f $(TARGETS) result_naive.txt result_strassen.txt
