#ifndef OPENMP_GS_PAGERANK_FUNCTIONS_H	/* Include guard */
#define OPENMP_GS_PAGERANK_FUNCTIONS_H

/* ===== INCLUDES ===== */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "coo_sparse_matrix.h"

/* ===== DEFINITIONS ===== */

//Colors used for better console output formating.
#define ANSI_COLOR_RED     "\x1B[31m"
#define ANSI_COLOR_GREEN   "\x1B[32m"
#define ANSI_COLOR_YELLOW  "\x1B[33m"
#define ANSI_COLOR_BLUE    "\x1B[34m"
#define ANSI_COLOR_CYAN    "\x1B[36m"
#define ANSI_COLOR_RESET   "\x1B[0m"

/* ===== CONSTANTS DEFINITION ===== */

// Constant strings that store the command line options available.
extern const char *ARGUMENT_CONVERGENCE_TOLERANCE;
extern const char *ARGUMENT_MAX_ITERATIONS;
extern const char *ARGUMENT_DAMPING_FACTOR;
extern const char *ARGUMENT_THREADS_NUMBER;
extern const char *ARGUMENT_VERBAL_OUTPUT;
extern const char *ARGUMENT_OUTPUT_HISTORY;
extern const char *ARGUMENT_OUTPUT_FILENAME;
// The numerical base used when parsing numerical command line arguments.
extern const int NUMERICAL_BASE;
// Default filename used for the output.
extern char *DEFAULT_OUTPUT_FILENAME;
// The size of the buffer used for reading the graph input file.
extern const int FILE_READ_BUFFER_SIZE;

/* ===== STRUCTURES ===== */

// A data structure to conveniently hold the algorithm's parameters.
typedef struct parameters {
	int numberOfPages, maxIterations;
	double convergenceCriterion, dampingFactor;
	bool verbose, history;
	char *outputFilename, *graphFilename;
} Parameters;

/* ===== FUNCTION DEFINITIONS ===== */

// Function validUsage outputs the correct way to use the program with command
// line arguments.
void validUsage(char *programName);

// Function checkIncrement is a helper function used in parseArguments (see
// bellow).
int checkIncrement(int previousIndex, int maxIndex, char *programName);

// Function parseArguments parses command line arguments.
void parseArguments(int argumentCount, char **argumentVector,
	Parameters *parameters);

// Function generateNormalizedTransitionMatrixFromFile reads through the entries
// of the file specified in the arguments (parameters->graphFilename), using
// them to populate the sparse array (transitionMatrix). The entries of the file
// represent the edges of the web transition graph. The entries are then
// modified to become the rows of the transition matrix.
void generateNormalizedTransitionMatrixFromFile(CsrSparseMatrix *transitionMatrix,
	Parameters *parameters);

// Function savePagerankToFile appends or overwrites the pagerank vector
// "pagerankVector" to the file with the filename supplied in the arguments.
void savePagerankToFile(char *filename, int *iterationsUntilConvergence,
	double *pagerankVector, int vectorSize, int iteration);

// Function initialize allocates memory for the pagerank vector, reads the
// dataset from the file and creates the transition probability distribution
// matrix.
void initialize(CsrSparseMatrix *transitionMatrix, double **pagerankVector,
	Parameters *parameters);

// Function vectorNorm calculates the first norm of a vector.
double vectorNorm(double *vector, int vectorSize);

// Function calculateNextPagerank calculates the next pagerank vector.
void calculateNextPagerank(CsrSparseMatrix *transitionMatrix,
	double *previousPagerankVector, double **pagerankVector,
	double *linksFromConvergedPagesPagerankVector,
	double *convergedPagerankVector, int vectorSize, double dampingFactor);

// Function pagerank iteratively calculates the pagerank of each page until
// either the convergence criterion is met or the maximum number of iterations
// is reached.
int* pagerank(CsrSparseMatrix *transitionMatrix, double **pagerankVector,
	bool *convergenceStatus, Parameters parameters, int* maxIterationsForConvergence);

#endif	// OPENMP_GS_PAGERANK_FUNCTIONS_H