#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>
#include <mpi.h> 
#include <iomanip>

#include "report_utils.hpp"

// --------------- read file (Chỉ dùng cho Master) -------------------
std::vector<double> readMatrix(const std::string& filename, int& rows, int& cols) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Can't open file " << filename << std::endl;
        exit(1);
    }
    infile >> rows >> cols;

    std::vector<double> matrix(rows * cols);
    for (long long i = 0; i < (long long)rows * cols; ++i) {
        infile >> matrix[i];
    }
    infile.close();
    return matrix;
}

// -------------- matrix multiplication -------------
void multiplyNaiveLocal(const std::vector<double>& localA, int localRows, int commonDim,
                        const std::vector<double>& B, int colsB,
                        std::vector<double>& localC) {
    
    // localRows: Số hàng A mà process này nắm giữ
    // commonDim: Số cột A (= Số hàng B)
    // colsB: Số cột B
    
    for (int i = 0; i < localRows; i++) {           
        for (int j = 0; j < colsB; j++) {    
            double sum = 0.0;
            for (int k = 0; k < commonDim; k++) {   
                sum += localA[i * commonDim + k] * B[k * colsB + j];
            }
            localC[i * colsB + j] = sum;
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
    // 1. Khoi tao MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <matrix1_file> <matrix2_file>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    std::string matrixFileA = argv[1];
    std::string matrixFileB = argv[2];

    int r1, c1, r2, c2;
    std::vector<double> A, B, C;

    // 2. Master (Rank 0) doc du lieu
    if (rank == 0) {
        A = readMatrix(matrixFileA, r1, c1);
        B = readMatrix(matrixFileB, r2, c2);

        if (c1 != r2) {
            std::cerr << "Size Error: Cols A != Rows B" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        std::cout << "Master: Doc xong (" << r1 << "x" << c1 << ") va (" << r2 << "x" << c2 << ")" << std::endl;
        C.resize(r1 * c2); // Cấp phát chỗ chứa kết quả cuối cùng
    }

    // 3. Broadcast kích thước cho tất cả các node con
    MPI_Bcast(&r1, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&c1, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r2, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&c2, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 4. Chuẩn bị dữ liệu cho Worker
    // Tất cả worker cần full ma trận B để nhân
    if (rank != 0) {
        B.resize(r2 * c2);
    }
    MPI_Bcast(B.data(), r2 * c2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // 5. Chia việc (Scattering A)
    // Tính số hàng mỗi process phải làm
    int rows_per_proc = r1 / size;
    int remainder = r1 % size;

    // Tính toán số hàng cụ thể cho process hiện tại (xử lý phần dư)
    int my_rows = (rank < remainder) ? rows_per_proc + 1 : rows_per_proc;
    
    // Cấp phát bộ nhớ cục bộ
    std::vector<double> localA(my_rows * c1);
    std::vector<double> localC(my_rows * c2);

    // Gửi dữ liệu từ Master sang Worker (Scatter thủ công để xử lý phần dư dễ hơn)
    if (rank == 0) {
        // Copy phần của Master
        std::copy(A.begin(), A.begin() + my_rows * c1, localA.begin());

        // Gửi cho các worker khác
        int offset = my_rows;
        for (int p = 1; p < size; p++) {
            int p_rows = (p < remainder) ? rows_per_proc + 1 : rows_per_proc;
            int count = p_rows * c1;
            MPI_Send(&A[offset * c1], count, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
            offset += p_rows;
        }
    } else {
        // Worker nhận phần của mình
        MPI_Recv(localA.data(), my_rows * c1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // 6. TÍNH TOÁN
    MPI_Barrier(MPI_COMM_WORLD); // Đồng bộ trước khi đo giờ
    auto start = std::chrono::high_resolution_clock::now();

    multiplyNaiveLocal(localA, my_rows, c1, B, c2, localC);

    MPI_Barrier(MPI_COMM_WORLD); // Đợi mọi người làm xong
    auto end = std::chrono::high_resolution_clock::now();

    // 7. Thu thập kết quả (Gathering C)
    if (rank == 0) {
        // Copy kết quả của Master vào C tổng
        std::copy(localC.begin(), localC.end(), C.begin());

        // Nhận từ các worker
        int offset = my_rows;
        for (int p = 1; p < size; p++) {
            int p_rows = (p < remainder) ? rows_per_proc + 1 : rows_per_proc;
            int count = p_rows * c2;
            MPI_Recv(&C[offset * c2], count, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            offset += p_rows;
        }

        // 8. Báo cáo (Chỉ Master làm)
        std::chrono::duration<double> elapsed = end - start;
        double checksum = calculateChecksum(C);

        std::cout << "Completion." << std::endl;
        std::cout << "MPI Processes: " << size << std::endl;
        std::cout << "Execution time: " << elapsed.count() << " second." << std::endl;
        std::cout << "Checksum: " << std::fixed << std::setprecision(5) << checksum << std::endl;

        appendReport("Naive", "MPI", r1, c1, r2, c2, elapsed.count(), checksum, size, -1);

    } else {
        // Worker gửi kết quả về Master
        MPI_Send(localC.data(), my_rows * c2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}