import sys
import random

def generate_matrix_file(filename, rows, cols):
    print(f"Generating {filename} with size {rows}x{cols}...")
    with open(filename, 'w') as f:
        # Ghi kích thước dòng đầu tiên
        f.write(f"{rows} {cols}\n")
        
        # Ghi dữ liệu
        for _ in range(rows):
            line = [str(random.uniform(0, 100)) for _ in range(cols)]
            f.write(" ".join(line) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python gen_matrix.py <r1> <c1> <r2> <c2>")
        sys.exit(1)

    r1, c1 = int(sys.argv[1]), int(sys.argv[2])
    r2, c2 = int(sys.argv[3]), int(sys.argv[4])

    if c1 != r2:
        print("Error: c1 must equal r2 for matrix multiplication.")
        sys.exit(1)

    generate_matrix_file("matrix1", r1, c1)
    generate_matrix_file("matrix2", r2, c2)
    print("Done.")