SHELL := /bin/bash

# ============================================
# COMMANDS

CC = gcc
RM = rm -f
CFLAGS=-O0 -g -I.
OBJ=serial_gs_pagerank.o serial_gs_pagerank_functions.o
DEPS=serial_gs_pagerank_functions.h

# ==========================================
# TARGETS

EXECUTABLES = pagerank

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
