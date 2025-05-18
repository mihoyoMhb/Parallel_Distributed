# Compiler and flags
MPICC = mpicc
CFLAGS = -O3 -Wall -g
LDFLAGS = -lm

# Source files
SRC = main.c data_distribution.c row_sort.c column_sort.c shear_sort.c
OBJ = $(SRC:.c=.o)
HEADERS = data_distribution.h row_sort.h column_sort.h shear_sort.h

# Target executable
TARGET = shear_sort

all: $(TARGET)

# Compile and link
$(TARGET): $(OBJ)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile object files
%.o: %.c $(HEADERS)
	$(MPICC) $(CFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJ) $(TARGET)

# Run a test with 16x16 matrix
test: $(TARGET)
	mpirun -np 4 ./$(TARGET) 16 12345

# Run a larger test
test-large: $(TARGET)
	mpirun -np 4 ./$(TARGET) 1024 12345

# Run with different process counts
test-scaling: $(TARGET)
	@echo "Running with 1 process:"
	mpirun -np 1 ./$(TARGET) 512 12345
	@echo "Running with 2 processes:"
	mpirun -np 2 ./$(TARGET) 512 12345
	@echo "Running with 4 processes:"
	mpirun -np 4 ./$(TARGET) 512 12345
	@echo "Running with 8 processes:"
	mpirun -np 8 ./$(TARGET) 512 12345

.PHONY: all clean test test-large test-scaling