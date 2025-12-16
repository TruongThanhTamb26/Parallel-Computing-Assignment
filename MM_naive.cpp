#include <iostream>
#include <vector>
#include <fstream> 
#include <string>
#include <chrono>   

#include "report_utils.hpp"

// --------------- read file -------------------
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

// -------------- matrix multiplication -------------
void multiplyNaive(const std::vector<double>& A, int rowsA, int colsA,
                   const std::vector<double>& B, int rowsB, int colsB,
                   std::vector<double>& C) {
    
    // Check column A == row B
    if (colsA != rowsB) {
        std::cerr << "Size Error: Column A (" << colsA << ") != Row B (" << rowsB << ")" << std::endl;
        exit(1);
    }

    int N = rowsA;
    int M = colsA; 
    int P = colsB;
    
    for (int i = 0; i < N; i++) {           
        for (int j = 0; j < P; j++) {    
            // DOT PRODUCT   
            double sum = 0.0;
            for (int k = 0; k < M; k++) {   
                sum += A[i * M + k] * B[k * P + j];
            }
            C[i * P + j] = sum;
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

int main() {
    int r1, c1, r2, c2;

    // Read file
    std::vector<double> A = readMatrix("matrix1", r1, c1);
    std::vector<double> B = readMatrix("matrix2", r2, c2);

    // Check Error
    if (c1 != r2) {
        std::cerr << "Size Error" << std::endl;
        return 1;
    }

    std::vector<double> C(r1 * c2);

    std::cout << "Start multiply matrix (" << r1 << "x" << c1 << ") and (" << r2 << "x" << c2 << ")..." << std::endl;

    // Time measurement
    auto start = std::chrono::high_resolution_clock::now();

    multiplyNaive(A, r1, c1, B, r2, c2, C);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    double checksum = calculateChecksum(C);

    // Result
    std::cout << "Completion." << std::endl;
    std::cout << "Execution time: " << elapsed.count() << " second." << std::endl;
    std::cout << "Checksum: " << checksum << std::endl;
    
    appendReport("Naive", "Serial", r1, c1, r2, c2, elapsed.count(), checksum, 1, 1);

    return 0;
}