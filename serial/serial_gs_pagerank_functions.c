#include "serial_gs_pagerank_functions.h"

const char *CONVERGENCE_ARGUMENT = "-c";
const char *MAX_ITERATIONS_ARGUMENT = "-m";
const char *DAMPING_FACTOR_ARGUMENT = "-a";
const char *VERBAL_OUTPUT_ARGUMENT = "-v";
const int NUMERICAL_BASE = 10;

void validUsage(char *programName) {
	printf("%s [-c convergence] [-m max_iterations] [-a alpha] [-v] <graph_file>\
		\n-c convergence\
		\n\tthe convergence criterion\
		\n-m max_iterations\
		\n\tmaximum number of iterations to perform\
		\n-a alpha\
		\n\tthe damping factor\
		\n-v enable verbal output\
		\n", programName);
	exit(EXIT_FAILURE);
}

int checkIncrement(int previousIndex, int maxIndex, char *programName) {
	if (previousIndex == maxIndex) {
		validUsage(programName);
		exit(EXIT_FAILURE);
	}
	return ++previousIndex;
}

void parseArguments(int argumentCount, char **argumentVector, Parameters *parameters) {
	if (argumentCount < 2 || argumentCount > 10) {
		validUsage(argumentVector[0]);
	}

	(*parameters).numberOfPages = 0;
	(*parameters).maxIterations = 0;
	(*parameters).convergenceCriterion = 1;
	(*parameters).dampingFactor = 0.85;
	(*parameters).verbose = false;

	char *endPointer;
	int argumentIndex = 1;

	while (argumentIndex < argumentCount) {
		if (!strcmp(argumentVector[argumentIndex], CONVERGENCE_ARGUMENT)) {
			argumentIndex = checkIncrement(argumentIndex, argumentCount, argumentVector[0]);

			double convergenceInput = strtod(argumentVector[argumentIndex], &endPointer);
			if (convergenceInput == 0) {
				printf("Invalid convergence argument\n");
				exit(EXIT_FAILURE);
			}
			(*parameters).convergenceCriterion = convergenceInput;
		} else if (!strcmp(argumentVector[argumentIndex], MAX_ITERATIONS_ARGUMENT)) {
			argumentIndex = checkIncrement(argumentIndex, argumentCount, argumentVector[0]);

			size_t iterationsInput = strtol(argumentVector[argumentIndex], &endPointer, NUMERICAL_BASE);
			if (iterationsInput == 0 && endPointer) {
				printf("Invalid iterations argument\n");
				exit(EXIT_FAILURE);
			}
			(*parameters).maxIterations = iterationsInput;
		} else if (!strcmp(argumentVector[argumentIndex], DAMPING_FACTOR_ARGUMENT)) {
			argumentIndex = checkIncrement(argumentIndex, argumentCount, argumentVector[0]);

			double alphaInput = strtod(argumentVector[argumentIndex], &endPointer);
			if ((alphaInput == 0 || alphaInput > 1) && endPointer) {
				printf("Invalid alpha argument\n");
				exit(EXIT_FAILURE);
			}
			(*parameters).dampingFactor = alphaInput;
		} else if (!strcmp(argumentVector[argumentIndex], VERBAL_OUTPUT_ARGUMENT)) {
			(*parameters).verbose = true;
		} else if (argumentIndex == argumentCount - 1) {
			(*parameters).graphFilename = argumentVector[argumentIndex];
		} else {
			validUsage(argumentVector[0]);
			exit(EXIT_FAILURE);
		}
		++argumentIndex;
	}
}

void readGraphFromFile(int ***directedWebGraph, Parameters *parameters) {
	FILE *graphFile;

	// Opens the file for reading
	graphFile = fopen((*parameters).graphFilename, "r+");
	if (!graphFile) {
		printf("Error opening file \n");
		exit(EXIT_FAILURE);
	}

	// Reads the dimensions of the (square) array from the file
	int readChar, numberOfLines=0;
	while((readChar = fgetc(graphFile))) {
		// Breaks if end of file
		if (readChar == EOF) break;
		// Otherwise, if the character is a break line, adds one to the count of lines
		if (readChar == '\n') {
			++numberOfLines;
		}
	}

	if ((*parameters).verbose) {
		printf("Line count of file is %d \n", numberOfLines);
	}

	// Each line of the file represents one page of the graph
	(*parameters).numberOfPages = numberOfLines;
	rewind(graphFile);

	// Allocates memory and loads values into directedWebGraph (matrix A)
	// Allocates memory for the rows
	(*directedWebGraph) = (int **) malloc((*parameters).numberOfPages * sizeof(int *));

	for (int i=0; i<(*parameters).numberOfPages; ++i) {
		// Allocates memory for the columns of this row
		(*directedWebGraph)[i] = (int *) malloc((*parameters).numberOfPages * sizeof(int));
		// Reads values from the file
		for (int j=0; j<(*parameters).numberOfPages; ++j) {
			if (!fscanf(graphFile, "%d ", &(*directedWebGraph)[i][j])) {
				break;
			}
			//printf("directedWebGraph[%d][%d] = %d", i , j, (*directedWebGraph)[i][j]);
		}
	}

	fclose(graphFile);
}

void generateNormalizedTransitionMatrix(double ***transitionMatrix,
	int **directedWebGraph, Parameters parameters) {
	// Allocates memory for the transitionMatrix rows
	(*transitionMatrix) = (double **) malloc(parameters.numberOfPages * sizeof(double *));

	for (int i=0; i<parameters.numberOfPages; ++i) {
		// Allocates memory for this row's columns
		(*transitionMatrix)[i] = (double *) malloc(parameters.numberOfPages * sizeof(double));

		int pageOutdegree = 0;
		//Calculates the outdegree of this page
		for (int j=0; j<parameters.numberOfPages; ++j) {
			pageOutdegree += directedWebGraph[i][j];
		}
		for (int j=0; j<parameters.numberOfPages; ++j) {
			if (pageOutdegree == 0) {
				// Introduces random jumps from dangling nodes (P' = P + D)
				// This makes sure that there are no pages with zero outdegree.
				(*transitionMatrix)[i][j] = 1. / parameters.numberOfPages;
			} else {
				(*transitionMatrix)[i][j] = 1. / pageOutdegree;
			}
		}
	}
}

void makeIrreducible(double ***transitionMatrix, Parameters parameters) {
	// Manipulates the values of transitionMatrix to make it irreducible. A
	// uniform probability (1/number_of_pages) and no personalization are used
	// here.

	// Introduces teleportation (P'' = cP' + (1 - c)E)
	for (int i=0; i<parameters.numberOfPages; ++i) {
		for (int j=0; j<parameters.numberOfPages; ++j) {
			(*transitionMatrix)[i][j] =
			parameters.dampingFactor *(*transitionMatrix)[i][j] +
			(1 - parameters.dampingFactor) / parameters.numberOfPages;
		}
	}
}

void transposeMatrix(double ***matrix, int rows, int columns) {
	// Transposes the matrix
	// Rows become columns and vice versa

	double **tempArray = (double **) malloc(rows * sizeof(double *));
	for (int i=0; i<rows; ++i) {
		tempArray[i] = malloc(columns * sizeof(double));

		for (int j=0; j<columns; ++j) {
			tempArray[i][j] = (*matrix)[j][i];
		}
	}

	//double **pointerToFreeMemoryLater = *matrix;
	matrix = &tempArray;
	/*for (int i=0; i<rows; ++i) {
		free(pointerToFreeMemoryLater[i]);
	}
	free(pointerToFreeMemoryLater);*/
}

void initialize(int ***directedWebGraph, double ***transitionMatrix,
	double **pagerankVector, Parameters *parameters) {

	if ((*parameters).verbose) {
		printf("----- Reading graph from file -----\n");
	}
	readGraphFromFile(directedWebGraph, parameters);

	if ((*parameters).verbose) {
		printf("\n----- Running with parameters -----\
			\nNumber of pages: %d", (*parameters).numberOfPages);
		if (!(*parameters).maxIterations) {
			printf("\nMaximum number of iterations: inf");
		} else {
			printf("\nMaximum number of iterations: %d", (*parameters).maxIterations);
		}
		printf("\nConvergence criterion: %f\
			\nDamping factor: %f\
			\nGraph filename: %s\n", (*parameters).convergenceCriterion,
			(*parameters).dampingFactor, (*parameters).graphFilename);
	}

	// Allocates memory for the pagerank vector
	(*pagerankVector) = (double *) malloc((*parameters).numberOfPages * sizeof(double));
	for (int i=0; i<(*parameters).numberOfPages; ++i) {
		(*pagerankVector)[i] = 1. / (*parameters).numberOfPages;
	}

	generateNormalizedTransitionMatrix(transitionMatrix, *directedWebGraph, *parameters);
	makeIrreducible(transitionMatrix, *parameters);
	transposeMatrix(transitionMatrix, (*parameters).numberOfPages, (*parameters).numberOfPages);
}

double vectorFirstNorm(double *vector, int vectorSize) {
	double norm = 0;

	for (int i=0; i<vectorSize; ++i) {
		norm += vector[i];
	}

	return norm;
}

void nextProbabilityDistribution(double ***transitionMatrix, double *previousPagerankVector,
	double **newPagerankVector, Parameters parameters) {

	transposeMatrix(transitionMatrix, parameters.numberOfPages, parameters.numberOfPages);
	for (int i=0; i<parameters.numberOfPages; ++i) {
		double sum = 0;

		for (int j=0; j<parameters.numberOfPages; ++j) {
			sum += (*transitionMatrix)[i][j] * previousPagerankVector[j];
		}
		(*newPagerankVector)[i] = parameters.dampingFactor * sum;
	}

	double normDifference = vectorFirstNorm(previousPagerankVector, parameters.numberOfPages) -
	vectorFirstNorm((*newPagerankVector), parameters.numberOfPages);

	for (int i=0; i<parameters.numberOfPages; ++i) {
		(*newPagerankVector)[i] += normDifference / parameters.numberOfPages;
	}

	transposeMatrix(transitionMatrix, parameters.numberOfPages, parameters.numberOfPages);
}

int pagerank(double ***transitionMatrix, double **pagerankVector, Parameters parameters) {
	int iterations = 0;
	double delta,
	*vectorDifference = (double *) malloc(parameters.numberOfPages * sizeof(double)),
	*previousPagerankVector = (double *) malloc(parameters.numberOfPages * sizeof(double));

	if (parameters.verbose) {
		printf("\n----- Starting iterations -----\n");
	}

	do {
		memcpy(previousPagerankVector, *pagerankVector, parameters.numberOfPages * sizeof(double));

		nextProbabilityDistribution(transitionMatrix, previousPagerankVector, pagerankVector, parameters);

		for (int i=0; i<parameters.numberOfPages; ++i) {
			vectorDifference[i] = (*pagerankVector)[i] - previousPagerankVector[i];
		}
		delta = vectorFirstNorm(vectorDifference, parameters.numberOfPages);

		++iterations;
		printf("Iteration %d: delta = %f\n", iterations, delta);
	} while (delta > parameters.convergenceCriterion &&
		(parameters.maxIterations != 0 || iterations < parameters.maxIterations));

	return iterations;
}