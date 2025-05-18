#!/bin/bash

# Path to the CSR matrix file
MATRIX_FILE="/home/mihoyohb/Datas/Data/28924x28924_2043492.csr"

# Number of processes to test
PROCESSES=(1 2 4 8 16)

# Number of repetitions for each test
REPETITIONS=5

# Output file to save results
OUTPUT_FILE="performance_results.txt"
CSV_FILE="performance_results.csv"

# Extract problem size from the file name
PROBLEM_SIZE=$(basename "$MATRIX_FILE" | cut -d'_' -f1)

# Create the CSV file header
echo "Problem Size,Processes,Average Time (seconds),Efficiency" > $CSV_FILE

echo "Performance Test Results" > $OUTPUT_FILE
echo "=========================" >> $OUTPUT_FILE
echo "Date: $(date)" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# Variable to store time for 1 process (needed for efficiency calculation)
time_one_process=0

for p in "${PROCESSES[@]}"; do
    echo "Testing with $p process(es)..." | tee -a $OUTPUT_FILE
    total_time=0
    
    for ((i=1; i<=$REPETITIONS; i++)); do
        echo "  Run $i of $REPETITIONS..." | tee -a $OUTPUT_FILE
        
        # Run the program and extract the time
        output=$(mpiexec -n $p ./mmv "$MATRIX_FILE")
        time_value=$(echo "$output" | head -n 1 | awk '{print $1}')
        
        echo "    Time: $time_value seconds" | tee -a $OUTPUT_FILE
        total_time=$(echo "$total_time + $time_value" | bc)
    done
    
    # Calculate average
    average=$(echo "scale=6; $total_time / $REPETITIONS" | bc)
    echo "  Average time with $p process(es): $average seconds" | tee -a $OUTPUT_FILE
    echo "" >> $OUTPUT_FILE
    
    # Store time for 1 process to calculate efficiency
    if [ "$p" -eq 1 ]; then
        time_one_process=$average
    fi
    
    # Calculate efficiency
    if [ "$p" -eq 1 ]; then
        efficiency=1
    else
        efficiency=$(echo "scale=6; $time_one_process / ($p * $average)" | bc)
    fi
    
    # Append to CSV
    echo "$PROBLEM_SIZE,$p,$average,$efficiency" >> $CSV_FILE
done

echo "Test completed. Results saved to $OUTPUT_FILE and $CSV_FILE"