import csv
import random
import argparse

def generate_matrix_csv(rows, cols, output_filename="matrix_data.csv"):
    """
    Generates a CSV file with matrix dimensions and random integer data.

    Args:
        rows (int): The number of rows in the matrix.
        cols (int): The number of columns in the matrix.
        output_filename (str): The name of the CSV file to create.
    """
    matrix = [[random.randint(0, 999) for _ in range(cols)] for _ in range(rows)]

    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write dimensions first
        writer.writerow([rows, cols])
        # Write matrix data
        writer.writerows(matrix)
    print(f"Generated '{output_filename}' with a {rows}x{cols} matrix.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a CSV file for matrix data.")
    parser.add_argument("--rows", type=int, default=6, help="Number of rows for the matrix.")
    parser.add_argument("--cols", type=int, default=6, help="Number of columns for the matrix.")
    parser.add_argument("--output", type=str, default="matrix_data.csv", help="Output CSV filename.")
    
    args = parser.parse_args()

    # Ensure N and M are divisible by N_proc and M_proc later in C
    # For now, just generate. Let's assume a 2x2 processor grid for initial testing,
    # so rows and cols should be even.
    if args.rows % 2 != 0:
        print(f"Warning: Number of rows ({args.rows}) is not even. Adjusting to {args.rows + 1} for easier division by 2 processes.")
        args.rows +=1
    if args.cols % 2 != 0:
        print(f"Warning: Number of columns ({args.cols}) is not even. Adjusting to {args.cols + 1} for easier division by 2 processes.")
        args.cols +=1

    generate_matrix_csv(args.rows, args.cols, args.output) 