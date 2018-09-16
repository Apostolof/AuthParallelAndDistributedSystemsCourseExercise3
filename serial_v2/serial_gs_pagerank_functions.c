#include "serial_gs_pagerank_functions.h"

const char *ARGUMENT_CONVERGENCE_TOLERANCE = "-c";
const char *ARGUMENT_MAX_ITERATIONS = "-m";
const char *ARGUMENT_DAMPING_FACTOR = "-a";
const char *ARGUMENT_VERBAL_OUTPUT = "-v";
const char *ARGUMENT_OUTPUT_HISTORY = "-h";
const char *ARGUMENT_OUTPUT_FILENAME = "-o";

const int NUMERICAL_BASE = 10;
char *DEFAULT_OUTPUT_FILENAME = "pagerank_output";

// ==================== PAGERANK ====================

int pagerank(int ***transitionMatrix, double **pagerankVector, Parameters parameters) {
	int iterations = 0;
	double delta,
	*vectorDifference = (double *) malloc(parameters.numberOfPages * sizeof(double)),
	*previousPagerankVector = (double *) malloc(parameters.numberOfPages * sizeof(double)),
	*convergedPagerankVector = (double *) malloc(parameters.numberOfPages * sizeof(double)),
	**linksFromConvergedPages = (double **) malloc(parameters.numberOfPages * sizeof(double *)),
	*linksFromConvergedPagesPagerankVector = (double *) malloc(parameters.numberOfPages * sizeof(double));
	bool *converganceMatrix = (bool *) malloc(parameters.numberOfPages * sizeof(bool));

	for (int i=0; i<parameters.numberOfPages; ++i) {
		convergedPagerankVector[i] = 0;
		converganceMatrix[i] = false;
		linksFromConvergedPagesPagerankVector[i] = 0;

		linksFromConvergedPages[i] = (double *) malloc(parameters.numberOfPages * sizeof(double));
		for (int j=0; j<parameters.numberOfPages; ++j) {
			linksFromConvergedPages[i][j] = 0;
		}
	}

	if (parameters.verbose) {
		printf("\n----- Starting iterations -----\n");
	}

	do {
		memcpy(previousPagerankVector, *pagerankVector, parameters.numberOfPages * sizeof(double));

		matrixVectorMultiplication(*transitionMatrix, previousPagerankVector,
			linksFromConvergedPagesPagerankVector, convergedPagerankVector,
			pagerankVector, parameters.numberOfPages, parameters.dampingFactor);

		if (parameters.history) {
			savePagerankToFile(parameters.outputFilename, iterations != 0,
				*pagerankVector, parameters.numberOfPages);
		}

		for (int i=0; i<parameters.numberOfPages; ++i) {
			vectorDifference[i] = (*pagerankVector)[i] - previousPagerankVector[i];
		}
		delta = vectorNorm(vectorDifference, parameters.numberOfPages);

		if (!iterations % 10) {
			for (int i=0; i<parameters.numberOfPages; ++i) {
				double temp = fabs((*pagerankVector)[i] - previousPagerankVector[i]) / fabs(previousPagerankVector[i]);
				if (temp < parameters.convergenceCriterion){
					converganceMatrix[i] = true;
					convergedPagerankVector[i] = (*pagerankVector)[i];
				}
			}

			for (int i=0; i<parameters.numberOfPages; ++i) {
				if (converganceMatrix[i] == true) {
					for (int j=0; j<parameters.numberOfPages; ++j){
						if (converganceMatrix[j] == false){
							linksFromConvergedPages[i][j] = (*transitionMatrix)[i][j];
						}
						// Zeros out CN and CC sub-matrices
						(*transitionMatrix)[i][j] = 0;
						// Zeros out NC sub-matrix
						(*transitionMatrix)[j][i] = 0;
					}

					double sum = 0;
					for (int j=0; j<parameters.numberOfPages; ++j) {
						sum += linksFromConvergedPages[i][j] * (*pagerankVector)[j];
					}
					linksFromConvergedPagesPagerankVector[i] = sum;
				}
			}
		}

		++iterations;
		printf("Iteration %d: delta = %f\n", iterations, delta);
	} while (delta > parameters.convergenceCriterion &&
		(parameters.maxIterations == 0 || iterations < parameters.maxIterations));

	if (!parameters.history) {
		savePagerankToFile(parameters.outputFilename, false, *pagerankVector,
			parameters.numberOfPages);
	}

	return iterations;
}

// ==================== INITIALIZATION ====================

/*
 * initialize allocates required memory for arrays, reads the web graph from the
 * from the file and creates the initial transition probability distribution
 * matrix.
*/
void initialize(int ***directedWebGraph, int ***transitionMatrix,
	double **pagerankVector, Parameters *parameters) {

	// Reads web graph from file
	if ((*parameters).verbose) {
		printf("----- Reading graph from file -----\n");
	}
	readGraphFromFile(directedWebGraph, parameters);

	// Outputs the algorithm parameters to the console
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
	double webUniformProbability = 1. / (*parameters).numberOfPages;
	for (int i=0; i<(*parameters).numberOfPages; ++i) {
		(*pagerankVector)[i] = webUniformProbability;
	}

	// Generates the initial transition matrix (matrix P).
	generateNormalizedTransitionMatrix(transitionMatrix, *directedWebGraph, *parameters);
	// Transposes the transition matrix (P^T).
	transposeMatrix(transitionMatrix, (*parameters).numberOfPages, (*parameters).numberOfPages);
}

/*
 * generateNormalizedTransitionMatrix generates the normalized transition matrix
 * from the graph data (matrix P').
*/
void generateNormalizedTransitionMatrix(int ***transitionMatrix,
	int **directedWebGraph, Parameters parameters) {
	// Allocates memory for the transitionMatrix rows
	(*transitionMatrix) = (int **) malloc(parameters.numberOfPages * sizeof(int *));

	for (int i=0; i<parameters.numberOfPages; ++i) {
		// Allocates memory for this row's columns
		(*transitionMatrix)[i] = (int *) malloc(parameters.numberOfPages * sizeof(int));
		/*
		// Calculates the outdegree of this page
		int pageOutdegree = 0;
		for (int j=0; j<parameters.numberOfPages; ++j) {
			pageOutdegree += directedWebGraph[i][j];
		}

		// Populates this row of the transition matrix
		if (pageOutdegree != 0) {
			// Calculates the uniform probability once.
			double pageUniformProbability = 1. / pageOutdegree;

			for (int j=0; j<parameters.numberOfPages; ++j) {
				if (directedWebGraph[i][j] == 1){
					(*transitionMatrix)[i][j] = pageUniformProbability;
				} else {
					(*transitionMatrix)[i][j] = 0;
				}
			}
		} else {
			for (int j=0; j<parameters.numberOfPages; ++j) {
				(*transitionMatrix)[i][j] = 0;
			}
		}
		*/
		
	}
	for (int i=0; i<parameters.numberOfPages; ++i){
		for (int j=0; i<parameters.numberOfPages; ++i){
			(*transitionMatrix)[i][j] = directedWebGraph[i][j];
		}
	}
}

// ==================== MATH UTILS ====================

/*
 * matrixVectorMultiplication calculates the product of the multiplication
 * between a matrix and the a vector in a cheap way.
*/
void matrixVectorMultiplication(int **transitionMatrix, double *previousPagerankVector,
	double *linksFromConvergedPagesPagerankVector, double *convergedPagerankVector,
	double **pagerankVector, int vectorSize, double dampingFactor) {
	double webUniformProbability = 1. / vectorSize;

	for (int i=0; i<vectorSize; ++i) {
		double sum1 = 0;
		double sum2 = 0;

		for (int j=0; j<i; ++j) {
			sum1 += transitionMatrix[i][j] * previousPagerankVector[j];
		}
		for (int j=i; j<vectorSize; ++j) {
			sum2 += transitionMatrix[i][j] * (*pagerankVector)[j];
		}
		(*pagerankVector)[i] = dampingFactor * sum1+dampingFactor * sum2+(1-dampingFactor);
	}

	//double normDifference = vectorNorm(previousPagerankVector, vectorSize) -
	//vectorNorm(*pagerankVector, vectorSize);

	for (int i=0; i<vectorSize; ++i) {
		//(*pagerankVector)[i] += normDifference * webUniformProbability +
		//linksFromConvergedPagesPagerankVector[i] + convergedPagerankVector[i];
		(*pagerankVector)[i] += linksFromConvergedPagesPagerankVector[i] + convergedPagerankVector[i];
	}
}

/*
 * vectorNorm calculates the first norm of a vector.
*/
double vectorNorm(double *vector, int vectorSize) {
	double norm = 0.;

	for (int i=0; i<vectorSize; ++i) {
		norm += fabs(vector[i]);
	}

	return norm;
}

/*
 * transposeMatrix transposes the matrix passed (by reference) in the arguments.
*/
void transposeMatrix(int ***matrix, int rows, int columns) {
	// Transposes the matrix
	// Rows become columns and vice versa

	int **tempArray = (int **) malloc(columns * sizeof(int *));
	for (int i=0; i<columns; ++i) {
		tempArray[i] = malloc(rows * sizeof(int));

		for (int j=0; j<rows; ++j) {
			tempArray[i][j] = (*matrix)[j][i];
		}
	}

	// TODO free memory

	//double **pointerToFreeMemoryLater = *matrix;
	*matrix = tempArray;
	/*for (int i=0; i<rows; ++i) {
		free(pointerToFreeMemoryLater[i]);
	}
	free(pointerToFreeMemoryLater);*/
}

// ==================== PROGRAM INPUT AND OUTPUT UTILS ====================

/*
 * parseArguments parses the command line arguments given by the user.
*/
void parseArguments(int argumentCount, char **argumentVector, Parameters *parameters) {
	if (argumentCount < 2 || argumentCount > 10) {
		validUsage(argumentVector[0]);
	}

	(*parameters).numberOfPages = 0;
	(*parameters).maxIterations = 0;
	(*parameters).convergenceCriterion = 1;
	(*parameters).dampingFactor = 0.85;
	(*parameters).verbose = false;
	(*parameters).history = false;
	(*parameters).outputFilename = DEFAULT_OUTPUT_FILENAME;

	char *endPointer;
	int argumentIndex = 1;

	while (argumentIndex < argumentCount) {
		if (!strcmp(argumentVector[argumentIndex], ARGUMENT_CONVERGENCE_TOLERANCE)) {
			argumentIndex = checkIncrement(argumentIndex, argumentCount, argumentVector[0]);

			double convergenceInput = strtod(argumentVector[argumentIndex], &endPointer);
			if (convergenceInput == 0) {
				printf("Invalid convergence argument\n");
				exit(EXIT_FAILURE);
			}
			(*parameters).convergenceCriterion = convergenceInput;
		} else if (!strcmp(argumentVector[argumentIndex], ARGUMENT_MAX_ITERATIONS)) {
			argumentIndex = checkIncrement(argumentIndex, argumentCount, argumentVector[0]);

			size_t iterationsInput = strtol(argumentVector[argumentIndex], &endPointer, NUMERICAL_BASE);
			if (iterationsInput == 0 && endPointer) {
				printf("Invalid iterations argument\n");
				exit(EXIT_FAILURE);
			}
			(*parameters).maxIterations = iterationsInput;
		} else if (!strcmp(argumentVector[argumentIndex], ARGUMENT_DAMPING_FACTOR)) {
			argumentIndex = checkIncrement(argumentIndex, argumentCount, argumentVector[0]);

			double alphaInput = strtod(argumentVector[argumentIndex], &endPointer);
			if ((alphaInput == 0 || alphaInput > 1) && endPointer) {
				printf("Invalid alpha argument\n");
				exit(EXIT_FAILURE);
			}
			(*parameters).dampingFactor = alphaInput;
		} else if (!strcmp(argumentVector[argumentIndex], ARGUMENT_VERBAL_OUTPUT)) {
			(*parameters).verbose = true;
		} else if (!strcmp(argumentVector[argumentIndex], ARGUMENT_OUTPUT_HISTORY)) {
			(*parameters).history = true;
		} else if (!strcmp(argumentVector[argumentIndex], ARGUMENT_OUTPUT_FILENAME)) {
			argumentIndex = checkIncrement(argumentIndex, argumentCount, argumentVector[0]);

			if (fopen(argumentVector[argumentIndex], "w") == NULL) {
				printf("Invalid output filename. Reverting to default.\n");
				continue;
			}
			(*parameters).outputFilename = argumentVector[argumentIndex];
		} else if (argumentIndex == argumentCount - 1) {
			(*parameters).graphFilename = argumentVector[argumentIndex];
		} else {
			validUsage(argumentVector[0]);
			exit(EXIT_FAILURE);
		}
		++argumentIndex;
	}
}

/*
 * readGraphFromFile loads the file supplied in the command line arguments to an
 * array (directedWebGraph) that represents the graph.
*/
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
		printf("Line count of file is %d \n", numberOfLines + 1);
	}

	// Each line of the file represents one page of the graph
	(*parameters).numberOfPages = numberOfLines + 1;
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
		}
	}

	fclose(graphFile);
}

/*
 * validUsage outputs a message to the console that informs the user of the
 * correct (valid) way to use the program.
*/
void validUsage(char *programName) {
	printf("%s [-c convergence_criterion] [-m max_iterations] [-a alpha] [-v] [-h] [-o output_filename] <graph_file>\
		\n-c convergence_criterion\
		\n\tthe convergence tolerance criterion\
		\n-m max_iterations\
		\n\tmaximum number of iterations to perform\
		\n-a alpha\
		\n\tthe damping factor\
		\n-v enable verbal output\
		\n-h enable history output to file\
		\n-o output_filename\
		\n\tfilename and path for the output\
		\n", programName);
	exit(EXIT_FAILURE);
}

/*
 * checkIncrement is a helper function for parseArguments function.
*/
int checkIncrement(int previousIndex, int maxIndex, char *programName) {
	if (previousIndex == maxIndex) {
		validUsage(programName);
		exit(EXIT_FAILURE);
	}
	return ++previousIndex;
}

void savePagerankToFile(char *filename, bool append, double *pagerankVector,
	int vectorSize) {
	FILE *outputFile;

	if (append) {
		outputFile = fopen(filename, "a");
	} else {
		outputFile = fopen(filename, "w");
	}

	if (outputFile == NULL) {
		printf("Error while opening the output file.\n");
		return;
	}

	for (int i=0; i<vectorSize; ++i) {
		fprintf(outputFile, "%f ", pagerankVector[i]);
	}
	fprintf(outputFile, "\n");

	fclose(outputFile);
}