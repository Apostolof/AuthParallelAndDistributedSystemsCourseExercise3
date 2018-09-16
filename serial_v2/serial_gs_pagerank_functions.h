#ifndef SERIAL_GS_PAGERANK_FUNCTIONS_H	/* Include guard */
#define SERIAL_GS_PAGERANK_FUNCTIONS_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * Constant strings that store the command line options available.
*/
extern const char *ARGUMENT_CONVERGENCE_TOLERANCE;
extern const char *ARGUMENT_MAX_ITERATIONS;
extern const char *ARGUMENT_DAMPING_FACTOR;
extern const char *ARGUMENT_VERBAL_OUTPUT;
extern const char *ARGUMENT_OUTPUT_HISTORY;
extern const char *ARGUMENT_OUTPUT_FILENAME;

// This is the numerical base used when parsing the numerical command line
// arguments.
extern const int NUMERICAL_BASE;
// Default filename used for the output.
extern char *DEFAULT_OUTPUT_FILENAME;

// Declares a data structure to conveniently hold the algorithm's parameters.
typedef struct parameters {
	int numberOfPages, maxIterations;
	double convergenceCriterion, dampingFactor;
	bool verbose, history;
	char *outputFilename, *graphFilename;
} Parameters;

// Function validUsage outputs the correct way to use the program with command
// line arguments.
void validUsage(char *programName);

// Function checkIncrement is a helper function used in parseArguments (see
// bellow).
int checkIncrement(int previousIndex, int maxIndex, char *programName);

// Function parseArguments parses command line arguments.
void parseArguments(int argumentCount, char **argumentVector, Parameters *parameters);

// Function readGraphFromFile loads adjacency matrix, that represents the web
// graph, stored in the file provided in the command line arguments to the array
// directedWebGraph.
void readGraphFromFile(int ***directedWebGraph, Parameters *parameters);

// Function savePagerankToFile appends or overwrites the pagerank vector
// "pagerankVector" to the file with the filename supplied in the arguments 
void savePagerankToFile(char *filename, bool append, double *pagerankVector,
	int vectorSize);

// Function generateNormalizedTransitionMatrix generates the normalized
// transition matrix from the web graph data.
void generateNormalizedTransitionMatrix(int ***transitionMatrix,
	int **directedWebGraph, Parameters parameters);

// Function transposeMatrix transposes a matrix.
void transposeMatrix(int ***matrix, int rows, int columns);

// Function initialize allocates required memory for arrays, reads the dataset
// from the file and creates the transition probability distribution matrix.
void initialize(
	int ***directedWebGraph, /*This is matrix G (web graph)*/
	int ***transitionMatrix, /*This is matrix A (transition probability distribution matrix)*/
	double **pagerankVector, /*This is the resulting pagerank vector*/
	Parameters *parameters
	);

// Function vectorNorm calculates the first norm of a vector.
double vectorNorm(double *vector, int vectorSize);

// Function matrixVectorMultiplication calculates the product of the
// multiplication between a matrix and the a vector.
void matrixVectorMultiplication(int **transitionMatrix, double *previousPagerankVector,
			double *linksFromConvergedPagesPagerankVector, double *convergedPagerankVector,
			double **pagerankVector, int vectorSize, double dampingFactor);

// Function pagerank iteratively calculates the pagerank of each page until
// either the convergence criterion is met or the maximum number of iterations
// is reached.
int pagerank(int ***transitionMatrix, double **pagerankVector, Parameters parameters);

#endif	// SERIAL_GS_PAGERANK_FUNCTIONS_H