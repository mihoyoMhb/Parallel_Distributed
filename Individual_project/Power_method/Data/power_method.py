#!/usr/bin/env python3
import numpy as np
import time
from scipy import sparse
import os
from scipy.sparse.linalg import eigs

def load_csr_matrix(file_path):
    """
    Load a CSR matrix from a custom CSR file format.
    
    The file format is:
    - First line: rows cols nnz
    - Second line: row_ptr array
    - Third line: col_ind array
    - Fourth line: values array
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        # Parse matrix dimensions and non-zero elements
        dimensions = lines[0].strip().split()
        rows, cols, nnz = map(int, dimensions)
        
        # Parse row_ptr, col_ind, and values arrays
        row_ptr = np.array(list(map(int, lines[1].strip().split())))
        col_ind = np.array(list(map(int, lines[2].strip().split())))
        data = np.array(list(map(float, lines[3].strip().split())))
        
        # Create CSR matrix
        matrix = sparse.csr_matrix((data, col_ind, row_ptr), shape=(rows, cols))
        
        return matrix, rows, cols, nnz

def main():
    # Get the current script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # File path for the CSR matrix
    csr_file = os.path.join(script_dir, "10974x10974_428650.csr")
    
    print(f"Loading CSR matrix from {csr_file}...")
    start_time = time.time()
    
    # Load the CSR matrix
    matrix, rows, cols, nnz = load_csr_matrix(csr_file)
    
    load_time = time.time() - start_time
    print(f"Matrix loaded in {load_time:.2f} seconds.")
    print(f"Matrix dimensions: {rows}x{cols} with {nnz} non-zero elements")
    
    # Run power method
    print("Running scipy.sparse.linalg.eigs to find dominant eigenvalue...")
    start_time = time.time()
    
    # Use eigs to find the largest eigenvalue (LM = Largest Magnitude)
    # eigs returns eigenvalues and eigenvectors. We ask for k=1 eigenvalue.
    # 'which='LM'' means largest magnitude.
    # We take the real part as eigenvalues can be complex.
    eigenvalues, eigenvectors = eigs(matrix, k=1, which='LM')
    eigenvalue = eigenvalues[0].real
    eigenvector = eigenvectors[:, 0].real # Eigenvectors are columns
    
    compute_time = time.time() - start_time
    
    print(f"Computation completed in {compute_time:.2f} seconds.")
    print(f"Dominant eigenvalue: {eigenvalue}")
    
    # Verify the result
    print("Verifying the result...")
    residual = np.linalg.norm(matrix.dot(eigenvector) - eigenvalue * eigenvector)
    print(f"Residual |Ax - λx|: {residual}")

if __name__ == "__main__":
    main() 