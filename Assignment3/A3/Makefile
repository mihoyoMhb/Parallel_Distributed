# Compiler
CC = mpicc


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


%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# "clean" target to remove generated files
clean:
	rm -f $(OBJS) $(TARGET) *.o core # Remove object files, target, any stray .o, and core dumps

# Phony targets are targets that don't represent actual files
.PHONY: all clean