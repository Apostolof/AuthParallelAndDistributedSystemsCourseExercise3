#ifndef SERIAL_GS_PAGERANK_FUNCTIONS_H	/* Include guard */
#define SERIAL_GS_PAGERANK_FUNCTIONS_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Constant strings that store the command line options available.
*/
extern const char *CONVERGENCE_ARGUMENT;
extern const char *MAX_ITERATIONS_ARGUMENT;
extern const char *DAMPING_FACTOR_ARGUMENT;
extern const char *VERBAL_OUTPUT_ARGUMENT;

// This is the numerical base used when parsing the numerical command line
// arguments.
extern const int NUMERICAL_BASE;

// Declares a data structure to conveniently hold the algorithm's parameters
typedef struct parameters {
	int numberOfPages, maxIterations;
	double convergenceCriterion, dampingFactor;
	bool verbose;
	char* graphFilename;
} Parameters;

// Function validUsage outputs the correct way to use the program with command
// line arguments.
void validUsage(char *programName);

// Function checkIncrement is a helper function used in parseArguments (see
// bellow).
int checkIncrement(int previousIndex, int maxIndex, char *programName);

// Function parseArguments parses command line arguments.
void parseArguments(int argumentCount, char **argumentVector, Parameters *parameters);

// Function readGraphFromFile loads the graph stored in the file provided in the
// command line arguments to the array directedWebGraph.
void readGraphFromFile(int ***directedWebGraph, Parameters *parameters);

// Function generateNormalizedTransitionMatrix generates the normalized transition
// matrix from the graph data.
void generateNormalizedTransitionMatrix(double ***transitionMatrix,
	int **directedWebGraph, Parameters parameters);

// Function makeIrreducible introduces teleportation to the transition matrix,
// making it irreducible.
void makeIrreducible(double ***transitionMatrix, Parameters parameters);

// Function transposeMatrix transposes a matrix.
void transposeMatrix(double ***matrix, int rows, int columns);

// Function initialize allocates required memory for arrays, reads the dataset
// from the file and creates the transition probability distribution matrix.
void initialize(
	int ***directedWebGraph, /*This is matrix G (web graph)*/
	double ***transitionMatrix, /*This is matrix A (transition probability distribution matrix)*/
	double **pagerankVector, /*This is the resulting pagerank vector*/
	Parameters *parameters
	);

// Function vectorFirstNorm calculates the first norm of a vector.
double vectorFirstNorm(double *vector, int vectorSize);

// Function nextProbabilityDistribution calculates the product of the transition
// matrix and the pagerank vector.
void nextProbabilityDistribution(double ***transitionMatrix, double *previousPagerankVector,
	double **newPagerankVector, Parameters parameters);

int pagerank(double ***transitionMatrix, double **pagerankVector, Parameters parameters);

#endif	// SERIAL_GS_PAGERANK_FUNCTIONS_H