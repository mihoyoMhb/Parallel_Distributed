# Compiler and flags
MPICC = mpicc
CFLAGS = -O3 -Wall -g -std=c99 # Added -std=c99 for compatibility like ceil in math.h
LDFLAGS = -lm

# Source files for V2
SRC_V2 = main_v2.c data_distribution_v2.c row_sort_v2.c column_sort_v2.c shear_sort_v2.c
OBJ_V2 = $(SRC_V2:.c=.o)
HEADERS_V2 = data_distribution_v2.h row_sort_v2.h column_sort_v2.h shear_sort_v2.h

# Target executable for V2
TARGET_V2 = shear_sort_v2

all: $(TARGET_V2)

# Compile and link V2
$(TARGET_V2): $(OBJ_V2)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile V2 object files
%.o: %.c $(HEADERS_V2)
	$(MPICC) $(CFLAGS) -c $< -o $@

# Target to run the V2 sort
# Assumes matrix_data.csv (or other specified CSV) exists in the current directory.
run: $(TARGET_V2)
	@echo "Running V2 Shear Sort with 4 processes..."
	@echo "Please ensure matrix_data.csv (or the CSV specified as an argument) exists."
	mpirun -np 4 ./$(TARGET_V2) matrix_data.csv # Defaulting to matrix_data.csv for run target

run_np2: $(TARGET_V2)
	@echo "Running V2 Shear Sort with 2 processes..."
	@echo "Please ensure matrix_data.csv (or the CSV specified as an argument) exists."
	mpirun -np 2 ./$(TARGET_V2) matrix_data.csv

run_np1: $(TARGET_V2)
	@echo "Running V2 Shear Sort with 1 process..."
	@echo "Please ensure matrix_data.csv (or the CSV specified as an argument) exists."
	mpirun -np 1 ./$(TARGET_V2) matrix_data.csv

# Clean up V2 files
clean:
	rm -f $(OBJ_V2) $(TARGET_V2)

.PHONY: all clean run run_np2 run_np1 