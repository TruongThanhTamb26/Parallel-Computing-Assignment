#include <iostream>
#include <vector>
#include <fstream> 
#include <string>
#include <chrono>   
#include <cmath>
#include <algorithm>

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

    // Divide size
    int newSize = size / 2;
    int subIdxSize = newSize * newSize;

    // fomulas
    /*
        M1 = (A11 + A22) (B11 + B22)
        M2 = (A21 + A22)B11
        M3 = A11(B12 − B22)
        M4 = A22(B21 − B11)
        M5 = (A11 + A12)B22
        M6 = (A21 − A11)(B11 + B12)
        M7 = (A12 − A22)(B21 + B22)

        C11 = M1 + M4 − M5 + M7
        C12 = M3 + M5
        C21 = M2 + M4
        C22 = M1 − M2 + M 3 + M6
    */

    vector<double> A11(subIdxSize), A12(subIdxSize), A21(subIdxSize), A22(subIdxSize);
    vector<double> B11(subIdxSize), B12(subIdxSize), B21(subIdxSize), B22(subIdxSize);
    vector<double> C11(subIdxSize), C12(subIdxSize), C21(subIdxSize), C22(subIdxSize);
    vector<double> M1(subIdxSize), M2(subIdxSize), M3(subIdxSize), M4(subIdxSize), M5(subIdxSize), M6(subIdxSize), M7(subIdxSize);
    vector<double> tempA(subIdxSize), tempB(subIdxSize);

    // Divide matrix
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

    // **Strassen algorithm
    
    // M1 = (A11 + A22) * (B11 + B22)
    add(A11, A22, tempA, newSize);
    add(B11, B22, tempB, newSize);
    strassenRecursive(tempA, tempB, M1, newSize);

    // M2 = (A21 + A22) * B11
    add(A21, A22, tempA, newSize);
    strassenRecursive(tempA, B11, M2, newSize);

    // M3 = A11 * (B12 - B22)
    sub(B12, B22, tempB, newSize);
    strassenRecursive(A11, tempB, M3, newSize);

    // M4 = A22 * (B21 - B11)
    sub(B21, B11, tempB, newSize);
    strassenRecursive(A22, tempB, M4, newSize);

    // M5 = (A11 + A12) * B22
    add(A11, A12, tempA, newSize);
    strassenRecursive(tempA, B22, M5, newSize);

    // M6 = (A21 - A11) * (B11 + B12)
    sub(A21, A11, tempA, newSize);
    add(B11, B12, tempB, newSize);
    strassenRecursive(tempA, tempB, M6, newSize);

    // M7 = (A12 - A22) * (B21 + B22)
    sub(A12, A22, tempA, newSize);
    add(B21, B22, tempB, newSize);
    strassenRecursive(tempA, tempB, M7, newSize);

    // 5. Tính C11, C12, C21, C22 từ M
    // C11 = M1 + M4 - M5 + M7
    add(M1, M4, tempA, newSize);
    sub(tempA, M5, tempB, newSize);
    add(tempB, M7, C11, newSize);

    // C12 = M3 + M5
    add(M3, M5, C12, newSize);

    // C21 = M2 + M4
    add(M2, M4, C21, newSize);

    // C22 = M1 - M2 + M3 + M6
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

// Wrapper function for handle Padding 
void multiplyStrassen(const vector<double>& A, int rowsA, int colsA,
                      const vector<double>& B, int rowsB, int colsB,
                      vector<double>& C) {
    
    // max size
    int maxSize = max({rowsA, colsA, rowsB, colsB});
    
    // Next power of 2
    int n = 1;
    while (n < maxSize) n *= 2;

    std::cout << "Strassen Padding: Resize from " << maxSize << " to " << n << "..." << std::endl;

    // Padded Matrices
    vector<double> A_pad(n * n, 0.0);
    vector<double> B_pad(n * n, 0.0);
    vector<double> C_pad(n * n, 0.0);

    // Copy data
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

    strassenRecursive(A_pad, B_pad, C_pad, n);

    // Unpadding
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

    std::cout << "Start multiply matrix with Strassen algorithm (" << r1 << "x" << c1 << ") and (" << r2 << "x" << c2 << ")..." << std::endl;

    // Time measurement
    auto start = std::chrono::high_resolution_clock::now();

    multiplyStrassen(A, r1, c1, B, r2, c2, C);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    double checksum = calculateChecksum(C);

    // Report
    std::cout << "Completion." << std::endl;
    std::cout << "Execution time: " << elapsed.count() << " second." << std::endl;
    std::cout << "Checksum: " << checksum << std::endl;

    std::ofstream outfile("result_strassen.txt");
    outfile << r1 << " " << c2 << "\n"; 
    for (int i = 0; i < r1; i++) {
        for (int j = 0; j < c2; j++) {
            outfile << C[i * c2 + j] << " ";
        }
        outfile << "\n";
    }
    outfile.close();

    return 0;
}