# Compiler and flags
CC = mpicc
CFLAGS = -std=c99 -O3 -Wall
LDFLAGS = -lm

# Target executable
TARGET = mmv

# Source files and object files
SRCS = main.c mmv.c
OBJS = $(SRCS:.c=.o)

# Default target
all: $(TARGET)

# Linking the executable
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# Compiling source files to object files
%.o: %.c mmv.h
	$(CC) $(CFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJS) $(TARGET)


.PHONY: all clean run run_with_matrix