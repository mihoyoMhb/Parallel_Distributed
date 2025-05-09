# Compiler
CC = mpicc

# Compiler flags
# -Wall: Enable all warnings
# -O2: Optimization level 2
# -g: Add debugging information (optional, useful for gdb)
# -std=c99 or -std=c11: Specify C standard (optional, good practice)
CFLAGS = -Wall -O3 -g -std=c99

# Linker flags (if any specific ones are needed, usually not for basic MPI)
LDFLAGS =

# Source files
SRCS = main.c quicksort.c pivot.c

# Object files (derived from source files)
OBJS = $(SRCS:.c=.o)

# Executable name
TARGET = quicksort

# Default target: build the executable
all: $(TARGET)

# Rule to link object files into the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Rule to compile .c files into .o files
# $< is the first prerequisite (the .c file)
# $@ is the target (the .o file)
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# "clean" target to remove generated files
clean:
	rm -f $(OBJS) $(TARGET) *.o core # Remove object files, target, any stray .o, and core dumps

# Phony targets are targets that don't represent actual files
.PHONY: all clean