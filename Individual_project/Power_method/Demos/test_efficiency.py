#!/usr/bin/env python3
import subprocess
import re
import os
import argparse
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import sys

# Configuration paths and names
C_SOURCE_FILE = "balanced_spmv.c"
EXECUTABLE_NAME = "balanced_spmv"

def list_csr_files():
    """List all CSR matrix files in the current directory"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    csr_files = [f for f in os.listdir(script_dir) if f.endswith('.csr')]
    
    if not csr_files:
        print("No CSR matrix files found. Please run mm_to_csr.py to convert matrix files first.")
        return None
    
    print("\nAvailable CSR matrix files:")
    for i, file in enumerate(csr_files, 1):
        # Extract matrix dimensions and non-zero elements from filename
        file_info = file.replace('.csr', '').split('_')
        if len(file_info) >= 2:
            dim = file_info[0]
            nnz = file_info[1]
            print(f"{i}. {file} (Size: {dim}, Non-zeros: {nnz})")
        else:
            print(f"{i}. {file}")
    
    return csr_files

def select_matrix_file():
    """Let the user select which matrix file to use"""
    csr_files = list_csr_files()
    
    if not csr_files:
        return None
    
    while True:
        try:
            choice = input("\nPlease select a matrix file number (enter q to quit): ")
            
            if choice.lower() == 'q':
                return None
            
            choice_idx = int(choice) - 1
            if 0 <= choice_idx < len(csr_files):
                selected_file = csr_files[choice_idx]
                print(f"Selected: {selected_file}")
                return selected_file
            else:
                print(f"Invalid selection. Please enter a number between 1 and {len(csr_files)}.")
        except ValueError:
            print("Please enter a valid number.")

def compile_c_code():
    """Compile C code using mpicc"""
    compile_command = [
        "mpicc", C_SOURCE_FILE,
        "-o", EXECUTABLE_NAME,
        "-lm",        # Link math library
        "-std=c99",   # Use C99 standard
        "-O3"         # Add optimization option
    ]
    print(f"Compiling C code: {' '.join(compile_command)}")
    
    try:
        process = subprocess.run(compile_command, check=True, capture_output=True, text=True)
        print("Compilation successful!")
        if process.stderr and process.stderr.strip():
            print(f"Compiler warnings:\n{process.stderr}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Compilation failed!")
        print(f"Return code: {e.returncode}")
        if e.stdout and e.stdout.strip():
            print(f"Standard output:\n{e.stdout}")
        if e.stderr and e.stderr.strip():
            print(f"Error output:\n{e.stderr}")
        return False
    except FileNotFoundError:
        print("Error: mpicc not found. Please make sure MPI is installed and in your PATH.")
        return False

def run_mpi_program(num_processes, matrix_file, max_retries=2):
    """Run the compiled MPI program and extract execution time"""
    executable_path = f"./{EXECUTABLE_NAME}"
    if not os.path.exists(executable_path):
        print(f"Error: Executable '{executable_path}' does not exist")
        return None
    
    # Run multiple times to get average, reducing random fluctuations
    times = []
    verification_status = "Unknown"
    
    for attempt in range(max_retries):
        command = ["mpiexec", "-n", str(num_processes), executable_path, matrix_file]
        print(f"Running MPI program (attempt {attempt+1}/{max_retries}): {' '.join(command)}")
        
        try:
            process = subprocess.run(command, check=True, capture_output=True, text=True, timeout=300)  # 5 minute timeout
            output = process.stdout
            
            # Extract computation time - UPDATED for English output
            time_match = re.search(r"Parallel computation time: (\d+\.\d+) seconds", output)
            if time_match:
                times.append(float(time_match.group(1)))
            
            # Check verification result - UPDATED for English output
            if "Verification passed" in output:
                verification_status = "Passed"
            elif "Verification failed" in output:
                verification_status = "Failed"
                mismatch_match = re.search(r"Verification failed with (\d+) mismatches", output)
                if mismatch_match:
                    verification_status += f" (Mismatches: {mismatch_match.group(1)})"
            
            # Only print complete output once
            if attempt == 0:
                print(f"Program output sample:\n{output[:500]}..." if len(output) > 500 else f"Program output:\n{output}")
        
        except subprocess.CalledProcessError as e:
            print(f"MPI program execution failed.")
            print(f"Return code: {e.returncode}")
            if e.stdout:
                print(f"Standard output:\n{e.stdout}")
            if e.stderr:
                print(f"Error output:\n{e.stderr}")
            return None
        except subprocess.TimeoutExpired:
            print(f"MPI program execution timed out (>5 minutes).")
            return None
        except FileNotFoundError:
            print(f"Error: mpiexec not found. Please make sure MPI is installed and in your PATH.")
            return None
    
    # Calculate average execution time
    if times:
        avg_time = sum(times) / len(times)
        print(f"Average execution time for {num_processes} processes: {avg_time:.6f} seconds")
        print(f"Verification result: {verification_status}")
        return avg_time
    else:
        print(f"Could not parse execution time for {num_processes} processes from output.")
        return None

def generate_plot(execution_times, baseline_p, output_file="speedup_efficiency.png"):
    """Generate speedup and efficiency plots"""
    if not execution_times:
        print("No execution time data available, cannot generate plots.")
        return
    
    baseline_time = execution_times.get(baseline_p)
    if not baseline_time:
        print(f"No baseline time (P={baseline_p}) available, cannot generate plots.")
        return
    
    sorted_p_counts = sorted(execution_times.keys())
    p_values = sorted_p_counts
    speedup_values = []
    efficiency_values = []
    
    for p in sorted_p_counts:
        time_p = execution_times[p]
        speedup = baseline_time / time_p if time_p > 0 else 0
        efficiency = (speedup / p) * baseline_p
        
        speedup_values.append(speedup)
        efficiency_values.append(efficiency)
    
    # Create plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Speedup plot
    ax1.plot(p_values, speedup_values, 'bo-', linewidth=2, markersize=8)
    ax1.plot(p_values, p_values, 'r--', linewidth=1, label='Ideal Speedup')
    ax1.set_title('Speedup vs. Process Count', fontsize=14)
    ax1.set_xlabel('Process Count (P)', fontsize=12)
    ax1.set_ylabel('Speedup (S)', fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    # Efficiency plot
    ax2.plot(p_values, efficiency_values, 'go-', linewidth=2, markersize=8)
    ax2.axhline(y=1.0, color='r', linestyle='--', label='Ideal Efficiency')
    ax2.set_title('Parallel Efficiency vs. Process Count', fontsize=14)
    ax2.set_xlabel('Process Count (P)', fontsize=12)
    ax2.set_ylabel('Parallel Efficiency (E)', fontsize=12)
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend()
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Test the performance and parallel efficiency of parallel sparse matrix-vector multiplication.")
    parser.add_argument("--processes", type=str, default="1,2,4,8", 
                        help="List of process counts to test (e.g., 1,2,4,8)")
    parser.add_argument("--no-plot", action="store_true", 
                        help="Do not generate performance plots")
    parser.add_argument("--retries", type=int, default=2,
                        help="Number of retries per process count to calculate average execution time")
    
    args = parser.parse_args()

    # Select matrix file to test
    matrix_file = select_matrix_file()
    if not matrix_file:
        print("No matrix file selected, test terminated.")
        return

    process_counts = []
    try:
        process_counts = [int(p.strip()) for p in args.processes.split(',') if p.strip()]
        if not process_counts or any(p <= 0 for p in process_counts):
            raise ValueError("Process counts must be positive integers.")
    except ValueError as e:
        print(f"Invalid --processes parameter: {e}. Example: --processes 1,2,4")
        return

    # Check if matrix file exists
    if not os.path.exists(matrix_file):
        print(f"Error: Matrix file '{matrix_file}' does not exist")
        return
    
    print("=" * 50)
    print(f"Starting Parallel SpMV Performance Test")
    print("=" * 50)
    print(f"Matrix file: {matrix_file}")
    print(f"Process counts to test: {process_counts}")
    print("=" * 50)

    # Compile C code
    if not compile_c_code():
        print("Test terminated due to compilation failure.")
        return

    execution_times = {}
    print("\n" + "=" * 50)
    print("Starting Performance Test Run")
    print("=" * 50)
    
    # Run tests for each process count
    for p_count in process_counts:
        print(f"\nTesting {p_count} processes...")
        time_taken = run_mpi_program(p_count, matrix_file, max_retries=args.retries)
        if time_taken is not None:
            execution_times[p_count] = time_taken
        else:
            print(f"Skipping efficiency calculation for {p_count} processes due to errors.")

    if not execution_times:
        print("No execution times recorded. Cannot calculate efficiency.")
        return

    # Determine baseline process count
    baseline_p = 1
    if 1 not in execution_times and execution_times:
        baseline_p = min(execution_times.keys())
        print(f"Baseline time for P=1 not available. Using P={baseline_p} as baseline.")
    elif 1 not in execution_times:
        print("Baseline time for P=1 not available and no other times recorded. Cannot calculate efficiency.")
        return
    
    baseline_time = execution_times.get(baseline_p)
    
    # Output efficiency table
    if baseline_time is not None:
        print("\n" + "=" * 50)
        print("Parallel Efficiency Results")
        print("=" * 50)
        print(f"Baseline time (P={baseline_p}): {baseline_time:.6f} seconds")
        print("-" * 60)
        print("| Process (P) | Time (sec) | Speedup (S) | Efficiency (E) |")
        print("|------------|------------|------------|--------------|")
        
        sorted_p_counts = sorted(execution_times.keys())

        for p_count in sorted_p_counts:
            time_p = execution_times[p_count]
            speedup = baseline_time / time_p if time_p > 0 else float('inf')
            efficiency = (speedup / p_count) * baseline_p
            print(f"| {p_count:<10} | {time_p:<10.4f} | {speedup:<10.2f} | {efficiency:<12.2f} |")
        
        print("-" * 60)
        
        # Generate performance plots
        if not args.no_plot:
            generate_plot(execution_times, baseline_p)
    else:
        print("Could not determine baseline time. Cannot print efficiency table.")

if __name__ == "__main__":
    main() 