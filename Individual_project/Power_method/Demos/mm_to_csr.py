#!/usr/bin/env python3
import sys
import os
import numpy as np
from scipy import sparse
import scipy.io as sio

def mm_to_csr(input_file, output_file):
    """Convert MatrixMarket format to custom CSR format"""
    print(f"Converting {input_file} to {output_file}...")
    
    try:
        # Read MatrixMarket format using scipy
        matrix = sio.mmread(input_file)
        
        # Convert to CSR format
        csr_matrix = sparse.csr_matrix(matrix)
        
        # Get matrix dimensions and number of non-zero elements
        rows, cols = csr_matrix.shape
        nnz = csr_matrix.nnz
        
        # Write to output file
        with open(output_file, 'w') as f:
            # First write dimensions and number of non-zero elements
            f.write(f"{rows} {cols} {nnz}\n")
            
            # Write row_ptr array (note CSR format uses 0-based indexing)
            for i in range(rows + 1):
                f.write(f"{csr_matrix.indptr[i]}")
                if i < rows:
                    f.write(" ")
            f.write("\n")
            
            # Write col_ind array
            for i in range(nnz):
                f.write(f"{csr_matrix.indices[i]}")
                if i < nnz - 1:
                    f.write(" ")
            f.write("\n")
            
            # Write values array
            for i in range(nnz):
                f.write(f"{csr_matrix.data[i]}")
                if i < nnz - 1:
                    f.write(" ")
            f.write("\n")
        
        print(f"Conversion complete! Matrix size: {rows}x{cols}, Non-zero elements: {nnz}")
        return rows, cols, nnz
    
    except Exception as e:
        print(f"Conversion failed: {str(e)}")
        return None

def process_all_mtx_files():
    """Process all mtx files in the Datas directory"""
    # Get script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Build Datas directory path
    datas_dir = os.path.join(script_dir, "Datas")
    
    if not os.path.exists(datas_dir):
        print(f"Error: Cannot find Datas directory ({datas_dir})")
        return False
    
    # Get all mtx files
    mtx_files = [f for f in os.listdir(datas_dir) if f.endswith('.mtx')]
    
    if not mtx_files:
        print(f"No .mtx files found in {datas_dir}")
        return False
    
    print(f"Found {len(mtx_files)} .mtx files, starting conversion...")
    success_count = 0
    
    for mtx_file in mtx_files:
        input_path = os.path.join(datas_dir, mtx_file)
        
        # First convert to get matrix information
        result = mm_to_csr(input_path, "temp.csr")
        
        if result:
            rows, cols, nnz = result
            # Name output file based on matrix dimensions and non-zero elements
            output_filename = f"{rows}x{cols}_{nnz}.csr"
            output_path = os.path.join(script_dir, output_filename)
            
            # Rename temporary file
            if os.path.exists("temp.csr"):
                os.rename("temp.csr", output_path)
                print(f"Created: {output_filename}")
                success_count += 1
    
    if success_count > 0:
        print(f"Successfully converted {success_count}/{len(mtx_files)} files")
        return True
    else:
        print("Failed to convert any files")
        return False

def main():
    if len(sys.argv) == 1:
        # No parameters, process all mtx files
        process_all_mtx_files()
    elif len(sys.argv) == 3:
        # Traditional usage to process a single file
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        
        if mm_to_csr(input_file, output_file):
            print(f"Successfully converted {input_file} to {output_file}")
        else:
            print("Conversion failed")
            sys.exit(1)
    else:
        print("Usage: python mm_to_csr.py [input_file.mtx output_file.csr]")
        print("If no parameters are provided, all .mtx files in the Datas directory will be processed")
        sys.exit(1)

if __name__ == "__main__":
    main() 