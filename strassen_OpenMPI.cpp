#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <string>
#include <vector>

#include "report_utils.hpp"

using namespace std;

struct TaskHeader {
    int blockSize; // dimension of the square sub-problem
    int taskId;    // identifier in {0..6}; -1 means stop
};

const int LEAF_SIZE = 64;
const int DISTRIBUTED_THRESHOLD = 256; // cutover for distributed Strassen

// --------------- READ FILE -------------------
std::vector<double> readMatrix(const std::string& filename, int& rows, int& cols) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Can't open file " << filename << std::endl;
        exit(1);
    }
    infile >> rows >> cols;
    std::cout << "Reading " << filename << " with " << rows << "x" << cols << "..." << std::endl;
    std::vector<double> matrix(rows * cols);
    for (long long i = 0; i < (long long)rows * cols; ++i) {
        infile >> matrix[i];
    }
    infile.close();
    return matrix;
}

// -------------- HELPER FUNCTIONS FOR STRASSEN -------------

// C = A + B
void add(const vector<double>& A, const vector<double>& B, vector<double>& C, int size) {
    for (int i = 0; i < size * size; i++) {
        C[i] = A[i] + B[i];
    }
}

// C = A - B
void sub(const vector<double>& A, const vector<double>& B, vector<double>& C, int size) {
    for (int i = 0; i < size * size; i++) {
        C[i] = A[i] - B[i];
    }
}

// C = A * B (naive)
void multiplyNaiveSquare(const vector<double>& A, const vector<double>& B, vector<double>& C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
}

void strassenRecursiveOMP(const vector<double>& A, const vector<double>& B, vector<double>& C, int size) {
    if (size <= LEAF_SIZE) {
        multiplyNaiveSquare(A, B, C, size);
        return;
    }

    int newSize = size / 2;
    int subIdxSize = newSize * newSize;

    vector<double> A11(subIdxSize), A12(subIdxSize), A21(subIdxSize), A22(subIdxSize);
    vector<double> B11(subIdxSize), B12(subIdxSize), B21(subIdxSize), B22(subIdxSize);
    vector<double> C11(subIdxSize), C12(subIdxSize), C21(subIdxSize), C22(subIdxSize);
    vector<double> M1(subIdxSize), M2(subIdxSize), M3(subIdxSize), M4(subIdxSize), M5(subIdxSize), M6(subIdxSize), M7(subIdxSize);
    vector<double> tempA(subIdxSize), tempB(subIdxSize);

    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            A11[i * newSize + j] = A[i * size + j];
            A12[i * newSize + j] = A[i * size + (j + newSize)];
            A21[i * newSize + j] = A[(i + newSize) * size + j];
            A22[i * newSize + j] = A[(i + newSize) * size + (j + newSize)];

            B11[i * newSize + j] = B[i * size + j];
            B12[i * newSize + j] = B[i * size + (j + newSize)];
            B21[i * newSize + j] = B[(i + newSize) * size + j];
            B22[i * newSize + j] = B[(i + newSize) * size + (j + newSize)];
        }
    }

    auto shouldSpawn = (newSize * newSize) >= 1024;

    #pragma omp task shared(M1) if(shouldSpawn)
    {
        vector<double> tA(subIdxSize), tB(subIdxSize);
        add(A11, A22, tA, newSize);
        add(B11, B22, tB, newSize);
        strassenRecursiveOMP(tA, tB, M1, newSize);
    }

    #pragma omp task shared(M2) if(shouldSpawn)
    {
        vector<double> tA(subIdxSize);
        add(A21, A22, tA, newSize);
        strassenRecursiveOMP(tA, B11, M2, newSize);
    }

    #pragma omp task shared(M3) if(shouldSpawn)
    {
        vector<double> tB(subIdxSize);
        sub(B12, B22, tB, newSize);
        strassenRecursiveOMP(A11, tB, M3, newSize);
    }

    #pragma omp task shared(M4) if(shouldSpawn)
    {
        vector<double> tB(subIdxSize);
        sub(B21, B11, tB, newSize);
        strassenRecursiveOMP(A22, tB, M4, newSize);
    }

    #pragma omp task shared(M5) if(shouldSpawn)
    {
        vector<double> tA(subIdxSize);
        add(A11, A12, tA, newSize);
        strassenRecursiveOMP(tA, B22, M5, newSize);
    }

    #pragma omp task shared(M6) if(shouldSpawn)
    {
        vector<double> tA(subIdxSize), tB(subIdxSize);
        sub(A21, A11, tA, newSize);
        add(B11, B12, tB, newSize);
        strassenRecursiveOMP(tA, tB, M6, newSize);
    }

    #pragma omp task shared(M7) if(shouldSpawn)
    {
        vector<double> tA(subIdxSize), tB(subIdxSize);
        sub(A12, A22, tA, newSize);
        add(B21, B22, tB, newSize);
        strassenRecursiveOMP(tA, tB, M7, newSize);
    }

    #pragma omp taskwait

    add(M1, M4, tempA, newSize);
    sub(tempA, M5, tempB, newSize);
    add(tempB, M7, C11, newSize);

    add(M3, M5, C12, newSize);

    add(M2, M4, C21, newSize);

    sub(M1, M2, tempA, newSize);
    add(tempA, M3, tempB, newSize);
    add(tempB, M6, C22, newSize);

    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            C[i * size + j] = C11[i * newSize + j];
            C[i * size + (j + newSize)] = C12[i * newSize + j];
            C[(i + newSize) * size + j] = C21[i * newSize + j];
            C[(i + newSize) * size + (j + newSize)] = C22[i * newSize + j];
        }
    }
}

void runStrassenOMP(const vector<double>& A, const vector<double>& B, vector<double>& C, int size) {
    if (size <= 0) return;
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            strassenRecursiveOMP(A, B, C, size);
        }
    }
}

void sendTask(int dest, int taskId, int blockSize, const vector<double>& left, const vector<double>& right) {
    TaskHeader header{blockSize, taskId};
    MPI_Send(&header, 2, MPI_INT, dest, 0, MPI_COMM_WORLD);
    int elems = blockSize * blockSize;
    MPI_Send(left.data(), elems, MPI_DOUBLE, dest, taskId, MPI_COMM_WORLD);
    MPI_Send(right.data(), elems, MPI_DOUBLE, dest, taskId, MPI_COMM_WORLD);
}

void sendStopSignal(int worldSize) {
    TaskHeader stop{-1, -1};
    for (int rank = 1; rank < worldSize; ++rank) {
        MPI_Send(&stop, 2, MPI_INT, rank, 0, MPI_COMM_WORLD);
    }
}

void workerLoop() {
    while (true) {
        TaskHeader header{};
        MPI_Status status;
        MPI_Recv(&header, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        if (header.blockSize <= 0 || header.taskId < 0) {
            break;
        }

        int elems = header.blockSize * header.blockSize;
        vector<double> left(elems), right(elems), result(elems, 0.0);
        MPI_Recv(left.data(), elems, MPI_DOUBLE, 0, header.taskId, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(right.data(), elems, MPI_DOUBLE, 0, header.taskId, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        runStrassenOMP(left, right, result, header.blockSize);
        MPI_Send(result.data(), elems, MPI_DOUBLE, 0, header.taskId, MPI_COMM_WORLD);
    }
}

void computeQuadrants(const vector<double>& src, vector<double>& q11, vector<double>& q12,
                      vector<double>& q21, vector<double>& q22, int size) {
    int newSize = size / 2;
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            q11[i * newSize + j] = src[i * size + j];
            q12[i * newSize + j] = src[i * size + (j + newSize)];
            q21[i * newSize + j] = src[(i + newSize) * size + j];
            q22[i * newSize + j] = src[(i + newSize) * size + (j + newSize)];
        }
    }
}

void combineBlocks(const vector<double>& C11, const vector<double>& C12,
                   const vector<double>& C21, const vector<double>& C22,
                   vector<double>& C, int size) {
    int newSize = size / 2;
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            C[i * size + j] = C11[i * newSize + j];
            C[i * size + (j + newSize)] = C12[i * newSize + j];
            C[(i + newSize) * size + j] = C21[i * newSize + j];
            C[(i + newSize) * size + (j + newSize)] = C22[i * newSize + j];
        }
    }
}

void distributedStrassen(const vector<double>& A_pad, const vector<double>& B_pad,
                          vector<double>& C_pad, int size, int worldSize) {
    if (worldSize == 1 || size <= DISTRIBUTED_THRESHOLD) {
        runStrassenOMP(A_pad, B_pad, C_pad, size);
        return;
    }

    int newSize = size / 2;
    int subIdxSize = newSize * newSize;

    vector<double> A11(subIdxSize), A12(subIdxSize), A21(subIdxSize), A22(subIdxSize);
    vector<double> B11(subIdxSize), B12(subIdxSize), B21(subIdxSize), B22(subIdxSize);

    computeQuadrants(A_pad, A11, A12, A21, A22, size);
    computeQuadrants(B_pad, B11, B12, B21, B22, size);

    vector<vector<double>> leftOps(7, vector<double>(subIdxSize));
    vector<vector<double>> rightOps(7, vector<double>(subIdxSize));
    vector<vector<double>> results(7, vector<double>(subIdxSize));

    add(A11, A22, leftOps[0], newSize);
    add(B11, B22, rightOps[0], newSize);

    add(A21, A22, leftOps[1], newSize);
    rightOps[1] = B11;

    leftOps[2] = A11;
    sub(B12, B22, rightOps[2], newSize);

    leftOps[3] = A22;
    sub(B21, B11, rightOps[3], newSize);

    add(A11, A12, leftOps[4], newSize);
    rightOps[4] = B22;

    sub(A21, A11, leftOps[5], newSize);
    add(B11, B12, rightOps[5], newSize);

    sub(A12, A22, leftOps[6], newSize);
    add(B21, B22, rightOps[6], newSize);

    for (int task = 0; task < 7; ++task) {
        int owner = task % worldSize;
        if (owner == 0) {
            runStrassenOMP(leftOps[task], rightOps[task], results[task], newSize);
        } else {
            sendTask(owner, task, newSize, leftOps[task], rightOps[task]);
            MPI_Recv(results[task].data(), subIdxSize, MPI_DOUBLE, owner, task, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    vector<double> tempA(subIdxSize), tempB(subIdxSize);
    vector<double> C11(subIdxSize), C12(subIdxSize), C21(subIdxSize), C22(subIdxSize);

    add(results[0], results[3], tempA, newSize);
    sub(tempA, results[4], tempB, newSize);
    add(tempB, results[6], C11, newSize);

    add(results[2], results[4], C12, newSize);
    add(results[1], results[3], C21, newSize);

    sub(results[0], results[1], tempA, newSize);
    add(tempA, results[2], tempB, newSize);
    add(tempB, results[5], C22, newSize);

    combineBlocks(C11, C12, C21, C22, C_pad, size);
}

void multiplyStrassenDistributed(const vector<double>& A, int rowsA, int colsA,
                                 const vector<double>& B, int rowsB, int colsB,
                                 vector<double>& C, int worldSize) {
    int maxSize = max({rowsA, colsA, rowsB, colsB});
    int n = 1;
    while (n < maxSize) n *= 2;

    std::cout << "Strassen Padding: Resize from " << maxSize << " to " << n << "..." << std::endl;

    vector<double> A_pad(n * n, 0.0);
    vector<double> B_pad(n * n, 0.0);
    vector<double> C_pad(n * n, 0.0);

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsA; j++) {
            A_pad[i * n + j] = A[i * colsA + j];
        }
    }
    for (int i = 0; i < rowsB; i++) {
        for (int j = 0; j < colsB; j++) {
            B_pad[i * n + j] = B[i * colsB + j];
        }
    }

    distributedStrassen(A_pad, B_pad, C_pad, n, worldSize);

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsB; j++) {
            C[i * colsB + j] = C_pad[i * n + j];
        }
    }
}

// ------------- check sum ---------------
double calculateChecksum(const std::vector<double>& matrix) {
    double sum = 0.0;
    for (double val : matrix) {
        sum += val;
    }
    return sum;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, worldSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    int threadCount = 1;
    #ifdef _OPENMP
        threadCount = omp_get_max_threads();
    #endif

    if (argc < 3) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <matrix1_file> <matrix2_file>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    if (rank != 0) {
        workerLoop();
        MPI_Finalize();
        return 0;
    }

    std::string matrixFileA = argv[1];
    std::string matrixFileB = argv[2];

    std::cout << "MPI Strassen. Processes: " << worldSize;
    #ifdef _OPENMP
        std::cout << " | OpenMP threads: " << threadCount;
    #endif
    std::cout << std::endl;

    int r1, c1, r2, c2;
    std::vector<double> A = readMatrix(matrixFileA, r1, c1);
    std::vector<double> B = readMatrix(matrixFileB, r2, c2);

    if (c1 != r2) {
        std::cerr << "Size Error" << std::endl;
        sendStopSignal(worldSize);
        MPI_Finalize();
        return 1;
    }

    std::vector<double> C(r1 * c2);

    std::cout << "Start distributed Strassen for matrices (" << r1 << "x" << c1
              << ") and (" << r2 << "x" << c2 << ")..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    multiplyStrassenDistributed(A, r1, c1, B, r2, c2, C, worldSize);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    double checksum = calculateChecksum(C);

    std::cout << "Completion." << std::endl;
    std::cout << "Execution time: " << elapsed.count() << " second." << std::endl;
    std::cout << "Checksum: " << checksum << std::endl;

    appendReport("Strassen", "Hybrid MPI+OpenMP", r1, c1, r2, c2, elapsed.count(), checksum, worldSize, threadCount);

    sendStopSignal(worldSize);
    MPI_Finalize();
    return 0;
}