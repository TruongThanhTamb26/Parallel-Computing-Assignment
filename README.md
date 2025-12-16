# Bài tập lớn: Nhân Ma trận Song song

Đây là bộ mã nguồn tham khảo cho bài tập lớn môn Lập trình Song song, bao gồm các cài đặt nhân ma trận đặc (dense matrix multiplication). Repository này chứa:

-   Nhân ma trận ngây thơ (Naive) $O(n^3)$ để kiểm tra tính đúng đắn và đo thời gian cơ sở.
-   Thuật toán Strassen chia để trị (phiên bản tuần tự, OpenMP, và MPI) với cơ chế tự động thêm lề (padding).
-   Các phiên bản nâng cấp sử dụng OpenMP và MPI để đánh giá khả năng mở rộng trên hệ thống chia sẻ bộ nhớ (shared-memory) và phân tán (distributed).
-   Script tiện ích để sinh ma trận đầu vào và ghi log kết quả.

---

## 📂 Cấu trúc Repository

| Đường dẫn | Mô tả |
| --- | --- |
| `MM_naive.cpp` | Cài đặt tuần tự 3 vòng lặp; ghi log chung vào `result_report.txt`. |
| `MM_strassen.cpp` | Strassen tuần tự với padding lũy thừa 2 và `LEAF_SIZE` tùy chỉnh; ghi log chung vào `result_report.txt`. |
| `naive_OpenMP.cpp`, `strassen_OpenMP.cpp` | Phiên bản chia sẻ bộ nhớ dùng OpenMP; báo cáo số luồng, thời gian và checksum vào `result_report.txt`. |
| `naive_MPI.cpp`, `strassen_MPI.cpp` | Phiên bản bộ nhớ phân tán dùng MPI; chỉ rank 0 ghi vào `result_report.txt`. |
| `naive_OpenMPI.cpp`, `strassen_OpenMPI.cpp` | Phiên bản Hybrid (MPI + OpenMP), dùng cùng file báo cáo `result_report.txt`. |
| `report_utils.hpp` | Hàm tiện ích `appendReport(...)` dùng chung để chuẩn hóa ghi log. |
| `gen_matrix.py` | Script sinh ma trận (`python3 gen_matrix.py r1 c1 r2 c2`). |
| `matrix1`, `matrix2` | File đầu vào mẫu định dạng text. |
| `Makefile` | File cấu hình biên dịch (`make`, `make openmp`, `make mpi`, `make gen`, `make clean`). |
| `result_report.txt` | File duy nhất chứa toàn bộ kết quả (% checksum, thời gian, số tiến trình/luồng). |

---

## ⚙️ Yêu cầu hệ thống (Prerequisites)

-   Linux hoặc môi trường POSIX tương tự (có `bash` và `make`).
-   Trình biên dịch C++17 (`g++`, `clang++`,...) có hỗ trợ OpenMP.
-   Môi trường MPI (`mpicxx`, `mpirun`).
-   Python 3 để sinh ma trận (khuyên dùng).

Cài đặt trên Ubuntu:

```bash
sudo apt install build-essential libomp-dev openmpi-bin openmpi-common python3
```

---

## 🧱 Biên dịch (Building)

| Mục tiêu | Lệnh | Ghi chú |
| --- | --- | --- |
| Tuần tự (Serial) | `make` | Tạo `MM_naive` và `MM_strassen`. |
| OpenMP | `make openmp` | Tạo `naive_OpenMP` và `strassen_OpenMP`. |
| MPI | `make mpi` | Tạo `naive_MPI` và `strassen_MPI`. |
| Tất cả | `make all openmp mpi` | Biên dịch toàn bộ. |

Các biến môi trường có thể thay đổi:

-   `THREADS`: Số luồng OpenMP dùng cho `make runopenmp` (mặc định `4`).
-   `PROCESS`: Số tiến trình MPI dùng cho `make runmpi` (mặc định `4`).

---

## ▶️ Hướng dẫn chạy (Running)

### Phiên bản Tuần tự

```bash
make run          # chạy MM_naive sau đó là MM_strassen với matrix1/matrix2
```

Mỗi chương trình sẽ in kích thước ma trận, thời gian thực thi, checksum, và thêm một block vào `result_report.txt` để dễ so sánh.

### Phiên bản OpenMP

```bash
THREADS=8 make runopenmp
```

Biến `OMP_NUM_THREADS` sẽ được tự động thiết lập. Bạn cũng có thể chạy thủ công:

```bash
export OMP_NUM_THREADS=8
./naive_OpenMP
./strassen_OpenMP
```

### Phiên bản MPI / Hybrid

```bash
PROCESS=6 make runmpi
```

Lệnh này tương đương với:

```bash
mpirun -np 6 ./naive_MPI
mpirun -np 6 ./strassen_MPI
```

Phiên bản MPI sẽ broadcast ma trận từ rank 0, chia việc cho các rank, thu thập kết quả và ghi log chung vào `result_report.txt`. Phiên bản `strassen_MPI` hiện tại phân tán mức đệ quy đầu tiên, trong khi bản `*_OpenMPI` tận dụng OpenMP cho tính toán cục bộ.

Sau khi chạy bất kỳ biến thể nào, mở `result_report.txt` để xem block tương ứng:

```
========================================
Algorithm : Naive
Mode      : MPI
Matrix A  : 1024 x 1024
Matrix B  : 1024 x 1024
Processes : 4
Threads   : 8          # chỉ xuất hiện khi có OpenMP
Time (s)  : 0.532871
Checksum  : 123456.000000

```

Các block được append theo thời gian, rất thuận tiện để so sánh tốc độ và checksum của các cấu hình khác nhau.

---

## 🗃️ Định dạng đầu vào & Sinh dữ liệu

Tất cả chương trình đều đọc ma trận dạng text:

```
<số hàng> <số cột>
v11 v12 ... v1c
...
vr1 vr2 ... vrc
```

Lệnh tạo dữ liệu mẫu:

```bash
# Tạo ma trận ngẫu nhiên 2048x2048 (giá trị 0..9)
python3 gen_matrix.py 2048 2048 2048 2048

# Hoặc dùng make (mặc định 500x500)
R1=1024 C1=1024 R2=1024 C2=1024 make gen
```

Đảm bảo kích thước hợp lệ (cột ma trận 1 == hàng ma trận 2).

---

## 🧮 Tổng quan Thuật toán

| Phiên bản | Đặc điểm |
| --- | --- |
| Naive (Serial/OpenMP/MPI) | 3 vòng lặp lồng nhau, không padding, dễ kiểm tra tính đúng đắn. |
| Strassen (Serial) | Padding lên lũy thừa 2, đệ quy đến `LEAF_SIZE` rồi chuyển sang nhân thường. |
| Strassen OpenMP | Tạo các task OpenMP cho 7 phép nhân con của Strassen và đồng bộ bằng `#pragma omp taskwait`. |
| Strassen MPI (Hybrid) | Rank 0 chia 7 bài toán con, gửi cho các rank khác. Các rank sử dụng OpenMP để tính toán song song cục bộ. |

Checksum (tổng tất cả phần tử) được in ra để đối chiếu kết quả giữa các phiên bản.

---

## 📊 Hướng dẫn Benchmark

1.  **Khởi động**: Chạy phiên bản tuần tự để lấy thời gian và checksum cơ sở.
2.  **Kiểm tra đúng đắn**: So sánh checksum của các bản song song với bản tuần tự.
3.  **Đánh giá mở rộng (Scaling)**: Thay đổi `THREADS` và `PROCESS`. Ghi lại bộ ba (kích thước, số tiến trình/luồng, thời gian).
4.  **Quét kích thước**: Thử nghiệm từ 100x100 đến 10.000x10.000. Strassen thường chỉ nhanh hơn Naive ở kích thước lớn.

---

## 🧹 Dọn dẹp

```bash
make clean              # xóa file thực thi
rm result_report.txt    # xóa toàn bộ báo cáo
```

---

## 🛠️ Xử lý sự cố (Troubleshooting)

-   **File not found**: Kiểm tra `matrix1`/`matrix2` có nằm cùng thư mục không.
-   **Dimension mismatch**: Kiểm tra dòng đầu tiên của file input.
-   **OpenMP không chạy song song**: Đảm bảo đã biên dịch với cờ `-fopenmp` và `OMP_NUM_THREADS` > 1.
-   **MPI bị treo**: Đảm bảo rank 0 luôn chạy đến cuối để gọi `MPI_Finalize`.
-   **Strassen chậm với ma trận nhỏ**: Giảm `LEAF_SIZE` hoặc dùng thuật toán Naive.

---

Chúc bạn làm bài tốt!