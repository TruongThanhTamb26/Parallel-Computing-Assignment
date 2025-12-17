import random

SIZES = [200, 1000, 4000]


def generate_matrix_file(filename, rows, cols):
    print(f"Generating {filename} with size {rows}x{cols}...")
    with open(filename, "w") as f:
        f.write(f"{rows} {cols}\n")
        for _ in range(rows):
            line = [str(random.uniform(0, 100)) for _ in range(cols)]
            f.write(" ".join(line) + "\n")


if __name__ == "__main__":
    for size in SIZES:
        filename1 = f"matrix1_{size}"
        filename2 = f"matrix2_{size}"
        generate_matrix_file(filename1, size, size)
        generate_matrix_file(filename2, size, size)

    print("Generated matrices for sizes:", ", ".join(map(str, SIZES)))