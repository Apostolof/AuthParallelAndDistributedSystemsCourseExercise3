SHELL := /bin/bash

# ============================================
# COMMANDS

CC = gcc
RM = rm -f
CFLAGS_DEBUG=-O0 -ggdb3 -Wall -I.
CFLAGS=-O3 -Wall -I.
OBJ=serial_gs_pagerank.o serial_gs_pagerank_functions.o coo_sparse_matrix.o csr_sparse_matrix.o
DEPS=serial_gs_pagerank_functions.h coo_sparse_matrix.h csr_sparse_matrix.h

# ==========================================
# TARGETS

EXECUTABLES = pagerank.out

.PHONY: all clean

all: $(EXECUTABLES)

# ==========================================
# DEPENDENCIES (HEADERS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PRECIOUS: $(EXECUTABLES) $(OBJ)

# ==========================================
# EXECUTABLE (MAIN)

$(EXECUTABLES): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	$(RM) *.o *~ $(EXECUTABLES)
