# Makefile
CC = mpicc
CFLAGS = -Wall -O3
TARGET = par_qsort
OBJS = main.o quicksort.o pivot.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean