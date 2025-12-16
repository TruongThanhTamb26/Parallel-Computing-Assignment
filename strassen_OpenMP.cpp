#include <iostream>
#include <vector>
#include <fstream> 
#include <string>
#include <chrono>   
#include <cmath>
#include <algorithm>
#include <omp.h> 

using namespace std;

const int LEAF_SIZE = 64; 

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


// -------------- STRASSEN RECURSIVE -------------
void strassenRecursive(const vector<double>& A, const vector<double>& B, vector<double>& C, int size) {
    // Base Case
    if (size <= LEAF_SIZE) {
        multiplyNaiveSquare(A, B, C, size);
        return;
    }

    int newSize = size / 2;
    int subIdxSize = newSize * newSize;

    // Allocate memory for sub-matrices
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
    
    #pragma omp task shared(A11, A22, B11, B22, M1) if(size > 128)
    {
        vector<double> tA(subIdxSize), tB(subIdxSize);
        add(A11, A22, tA, newSize);
        add(B11, B22, tB, newSize);
        strassenRecursive(tA, tB, M1, newSize);
    }

    #pragma omp task shared(A21, A22, B11, M2) if(size > 128)
    {
        vector<double> tA(subIdxSize);
        add(A21, A22, tA, newSize);
        strassenRecursive(tA, B11, M2, newSize);
    }

    #pragma omp task shared(A11, B12, B22, M3) if(size > 128)
    {
        vector<double> tB(subIdxSize);
        sub(B12, B22, tB, newSize);
        strassenRecursive(A11, tB, M3, newSize);
    }

    #pragma omp task shared(A22, B21, B11, M4) if(size > 128)
    {
        vector<double> tB(subIdxSize);
        sub(B21, B11, tB, newSize);
        strassenRecursive(A22, tB, M4, newSize);
    }

    #pragma omp task shared(A11, A12, B22, M5) if(size > 128)
    {
        vector<double> tA(subIdxSize);
        add(A11, A12, tA, newSize);
        strassenRecursive(tA, B22, M5, newSize);
    }

    #pragma omp task shared(A21, A11, B11, B12, M6) if(size > 128)
    {
        vector<double> tA(subIdxSize), tB(subIdxSize);
        sub(A21, A11, tA, newSize);
        add(B11, B12, tB, newSize);
        strassenRecursive(tA, tB, M6, newSize);
    }

    #pragma omp task shared(A12, A22, B21, B22, M7) if(size > 128)
    {
        vector<double> tA(subIdxSize), tB(subIdxSize);
        sub(A12, A22, tA, newSize);
        add(B21, B22, tB, newSize);
        strassenRecursive(tA, tB, M7, newSize);
    }

    // Wait for all tasks to complete before merging results
    #pragma omp taskwait

    add(M1, M4, tempA, newSize);
    sub(tempA, M5, tempB, newSize);
    add(tempB, M7, C11, newSize);

    add(M3, M5, C12, newSize);
    add(M2, M4, C21, newSize);

    sub(M1, M2, tempA, newSize);
    add(tempA, M3, tempB, newSize);
    add(tempB, M6, C22, newSize);

    // Merge
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            C[i * size + j] = C11[i * newSize + j];
            C[i * size + (j + newSize)] = C12[i * newSize + j];
            C[(i + newSize) * size + j] = C21[i * newSize + j];
            C[(i + newSize) * size + (j + newSize)] = C22[i * newSize + j];
        }
    }
}

// Wrapper function 
void multiplyStrassen(const vector<double>& A, int rowsA, int colsA,
                      const vector<double>& B, int rowsB, int colsB,
                      vector<double>& C) {
    
    int maxSize = max({rowsA, colsA, rowsB, colsB});
    int n = 1;
    while (n < maxSize) n *= 2;

    std::cout << "Strassen Padding: Resize from " << maxSize << " to " << n << "..." << std::endl;

    vector<double> A_pad(n * n, 0.0);
    vector<double> B_pad(n * n, 0.0);
    vector<double> C_pad(n * n, 0.0);

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsA; j++) A_pad[i * n + j] = A[i * colsA + j];
    }
    for (int i = 0; i < rowsB; i++) {
        for (int j = 0; j < colsB; j++) B_pad[i * n + j] = B[i * colsB + j];
    }

    #pragma omp parallel
    {
        #pragma omp single
        {
            strassenRecursive(A_pad, B_pad, C_pad, n);
        }
    }

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsB; j++) C[i * colsB + j] = C_pad[i * n + j];
    }
}

// ------------- check sum ---------------
double calculateChecksum(const std::vector<double>& matrix) {
    double sum = 0.0;
    for (double val : matrix) sum += val;
    return sum;
}

int main() {
    int r1, c1, r2, c2;
    
    #ifdef _OPENMP
        std::cout << "OpenMP Strassen. Max threads: " << omp_get_max_threads() << std::endl;
    #endif

    std::vector<double> A = readMatrix("matrix1", r1, c1);
    std::vector<double> B = readMatrix("matrix2", r2, c2);

    if (c1 != r2) {
        std::cerr << "Size Error" << std::endl;
        return 1;
    }

    std::vector<double> C(r1 * c2);
    std::cout << "Start multiply..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    multiplyStrassen(A, r1, c1, B, r2, c2, C);
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    double checksum = calculateChecksum(C);

    std::cout << "Completion." << std::endl;
    std::cout << "Execution time: " << elapsed.count() << " second." << std::endl;
    std::cout << "Checksum: " << checksum << std::endl;

    std::ofstream outfile("result_openmp.txt", std::ios::app);
    outfile << "Checksum of strassen: " << checksum << std::endl;
    outfile << "Execution time of strassen: " << elapsed.count() << " second." << std::endl;
    outfile.close();
    return 0;
}