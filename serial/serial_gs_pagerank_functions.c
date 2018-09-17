#include "serial_gs_pagerank_functions.h"

const char *ARGUMENT_CONVERGENCE_TOLERANCE = "-c";
const char *ARGUMENT_MAX_ITERATIONS = "-m";
const char *ARGUMENT_DAMPING_FACTOR = "-a";
const char *ARGUMENT_VERBAL_OUTPUT = "-v";
const char *ARGUMENT_OUTPUT_HISTORY = "-h";
const char *ARGUMENT_OUTPUT_FILENAME = "-o";

const int NUMERICAL_BASE = 10;
char *DEFAULT_OUTPUT_FILENAME = "pagerank_output";
const int MAX_PAGE_LINKS_TEXT_SIZE = 4096;

// ==================== PAGERANK ====================

int pagerank(SparseMatrix *transitionMatrix, double **pagerankVector, Parameters parameters) {
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

		matrixVectorMultiplication(transitionMatrix, previousPagerankVector,
			pagerankVector, parameters.numberOfPages, parameters.dampingFactor);

		for (int i=0; i<parameters.numberOfPages; ++i) {
			(*pagerankVector)[i] += linksFromConvergedPagesPagerankVector[i] + convergedPagerankVector[i];
		}

		if (parameters.history) {
			savePagerankToFile(parameters.outputFilename, iterations != 0,
				*pagerankVector, parameters.numberOfPages);
		}

		for (int i=0; i<parameters.numberOfPages; ++i) {
			vectorDifference[i] = (*pagerankVector)[i] - previousPagerankVector[i];
		}
		delta = vectorNorm(vectorDifference, parameters.numberOfPages);

		if (iterations && !iterations % 10) {
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
							SparseMatrixElement *element = getElement(*transitionMatrix, i, j);
							linksFromConvergedPages[i][j] = element != NULL ? element->value : 0;
						}
						deleteElement(transitionMatrix, i, j);
						deleteElement(transitionMatrix, j, i);
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
void initialize(SparseMatrix *transitionMatrix,
	double **pagerankVector, Parameters *parameters) {

	// Reads web graph from file
	if ((*parameters).verbose) {
		printf("----- Reading graph from file -----\n");
	}
	generateNormalizedTransitionMatrixFromFile(transitionMatrix, parameters);

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

	// Transposes the transition matrix (P^T).
	transposeSparseMatrix(transitionMatrix);
}

// ==================== MATH UTILS ====================

/*
 * matrixVectorMultiplication calculates the product of the multiplication
 * between a matrix and the a vector in a cheap way.
*/
void matrixVectorMultiplication(SparseMatrix *transitionMatrix, double *previousPagerankVector,
	double **pagerankVector, int vectorSize, double dampingFactor) {
	double webUniformProbability = 1. / vectorSize;

	sparseMatrixVectorMultiplication(*transitionMatrix, previousPagerankVector,
		pagerankVector, vectorSize);

	for (int i=0; i<vectorSize; ++i) {
		(*pagerankVector)[i] = dampingFactor * (*pagerankVector)[i];
	}

	double normDifference = vectorNorm(previousPagerankVector, vectorSize) -
	vectorNorm(*pagerankVector, vectorSize);

	for (int i=0; i<vectorSize; ++i) {
		(*pagerankVector)[i] += normDifference * webUniformProbability;
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
void generateNormalizedTransitionMatrixFromFile(SparseMatrix *transitionMatrix,
	Parameters *parameters){
	FILE *graphFile;

	// Opens the file for reading
	graphFile = fopen((*parameters).graphFilename, "r+");
	if (!graphFile) {
		printf("Error opening file \n");
		exit(EXIT_FAILURE);
	}

	int pageIndex, count = 0;
	while (fscanf(graphFile, "%d:", &pageIndex) != EOF) {
		if (!(pageIndex%51050)) {
			printf("\t%d\t%d%%\n", pageIndex, ++count);
		}

		char *restOfLine = malloc(MAX_PAGE_LINKS_TEXT_SIZE);
		if (!fgets(restOfLine, MAX_PAGE_LINKS_TEXT_SIZE, graphFile)) {
			exit(EXIT_FAILURE);
		}

		char *token = strtok(restOfLine, " ");

		while (token != NULL) {
			if (strcmp(token, "\n") == 0) {
				//token = strtok (NULL, " ");
				break;
			}

			int outLink = atoi(token);
			if (outLink != -1) {
				apendElement(transitionMatrix, 1, pageIndex, outLink);
			}
			token = strtok (NULL, " ");
		}
	}
	printf("\t100%%\n");
	printf("number of edges = %d\n", transitionMatrix->elements);

	(*parameters).numberOfPages = pageIndex + 1;

	int currentRow = transitionMatrix->firstElement->rowIndex;
	SparseMatrixElement *startElement = transitionMatrix->firstElement;
	while(true) {
		int pageOutdegree = 1;
		SparseMatrixElement *currentElement = startElement->nextElement;

		// Calculates current page's outdegree
		while (currentElement != NULL) {
			if (currentElement->rowIndex == currentRow) {
				++pageOutdegree;
				currentElement = currentElement->nextElement;
			} else {
				break;
			}
		}

		// Assigns the value 1/outdegree to current page's columns
		currentElement = startElement;
		for (int i=0; i<pageOutdegree; ++i) {
			if (currentElement->rowIndex == currentRow) {
				currentElement->value = 1. / pageOutdegree;
				currentElement = currentElement->nextElement;
			} else {
				break;
			}
		}

		// Reached the last element;
		if (currentElement == NULL) {
			break;
		}

		startElement = currentElement;
		currentRow = startElement->rowIndex;
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