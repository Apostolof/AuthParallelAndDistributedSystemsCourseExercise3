#include "serial_gs_pagerank_functions.h"

const char *ARGUMENT_CONVERGENCE_TOLERANCE = "-c";
const char *ARGUMENT_MAX_ITERATIONS = "-m";
const char *ARGUMENT_DAMPING_FACTOR = "-a";
const char *ARGUMENT_VERBAL_OUTPUT = "-v";
const char *ARGUMENT_OUTPUT_HISTORY = "-h";
const char *ARGUMENT_OUTPUT_FILENAME = "-o";

const int NUMERICAL_BASE = 10;
char *DEFAULT_OUTPUT_FILENAME = "pagerank_output";
const int FILE_READ_BUFFER_SIZE = 4096;

// ==================== PAGERANK ====================

int pagerank(SparseMatrix *transitionMatrix, double **pagerankVector,
	bool *convergenceStatus, Parameters parameters) {
	int iterations = 0;
	double delta,
	*vectorDifference = (double *) malloc(parameters.numberOfPages * sizeof(double)),
	*previousPagerankVector = (double *) malloc(parameters.numberOfPages * sizeof(double)),
	*convergedPagerankVector = (double *) malloc(parameters.numberOfPages * sizeof(double)),
	**linksFromConvergedPages = (double **) malloc(parameters.numberOfPages * sizeof(double *)),
	*linksFromConvergedPagesPagerankVector = (double *) malloc(parameters.numberOfPages * sizeof(double));
	bool *converganceMatrix = (bool *) malloc(parameters.numberOfPages * sizeof(bool));
	*convergenceStatus = false;

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
		printf(ANSI_COLOR_YELLOW "\n----- Starting iterations -----\n" ANSI_COLOR_RESET);
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
		if (delta < parameters.convergenceCriterion) {
			*convergenceStatus = true;
		}

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
		if (iterations%2) {
			printf(ANSI_COLOR_BLUE "Iteration %d: delta = %f\n" ANSI_COLOR_RESET, iterations, delta);
		} else {
			printf(ANSI_COLOR_CYAN "Iteration %d: delta = %f\n" ANSI_COLOR_RESET, iterations, delta);
		}
	} while (!*convergenceStatus &&
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
		printf(ANSI_COLOR_YELLOW "----- Reading graph from file -----\n" ANSI_COLOR_RESET);
	}
	generateNormalizedTransitionMatrixFromFile(transitionMatrix, parameters);

	// Outputs the algorithm parameters to the console
	if ((*parameters).verbose) {
		printf(ANSI_COLOR_YELLOW "\n----- Running with parameters -----\n" ANSI_COLOR_RESET\
			"Number of pages: %d", (*parameters).numberOfPages);
		if (!(*parameters).maxIterations) {
			printf("\nMaximum number of iterations: inf");
		} else {
			printf("\nMaximum number of iterations: %d", (*parameters).maxIterations);
		}
		printf("\nConvergence criterion: %f" \
			"\nDamping factor: %f" \
			"\nGraph filename: %s\n", (*parameters).convergenceCriterion,
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

	char buffer[FILE_READ_BUFFER_SIZE];
	char *readResult;
	// Skips the first two lines
	readResult = fgets(buffer, FILE_READ_BUFFER_SIZE, graphFile);
	readResult = fgets(buffer, FILE_READ_BUFFER_SIZE, graphFile);
	if (readResult == NULL) {
		printf("Error while reading from the file. Does the file have the correct format?\n");
		exit(EXIT_FAILURE);
	}

	// Third line contains the numbers of nodes and edges
	int numberOfNodes = 0, numberOfEdges;

	readResult = fgets(buffer, FILE_READ_BUFFER_SIZE, graphFile);
	if (readResult == NULL) {
		printf("Error while reading from the file. Does the file have the correct format?\n");
		exit(EXIT_FAILURE);
	}

	// Parses the number of nodes and number of edges
	{
		// Splits string to whitespace
		char *token = strtok(buffer, " ");
		bool nextIsNodes = false, nextIsEdges = false;

		while (token != NULL) {
			if (strcmp(token, "Nodes:") == 0) {
				nextIsNodes = true;
			} else if (nextIsNodes) {
				numberOfNodes = atoi(token);
				nextIsNodes = false;
			} else if (strcmp(token, "Edges:") == 0) {
				nextIsEdges = true;
			} else if (nextIsEdges) {
				numberOfEdges = atoi(token);
				break;
			}

			// Gets next string token
			token = strtok (NULL, " ,.-");
		}
	}

	if ((*parameters).verbose) {
		printf("The number of pages is: %d\nThe number of edges is: %d\n",
			numberOfNodes, numberOfEdges);
	}
	(*parameters).numberOfPages = numberOfNodes;

	// Skips the fourth line
	readResult = fgets(buffer, 512, graphFile);
	if (readResult == NULL) {
		printf("Error while reading from the file. Does the file have the correct format?\n");
		exit(EXIT_FAILURE);
	}

	printf("SIZE OF STRUCT = %lu Bytes\n", sizeof(SparseMatrixElement));

	int fivePercentIncrements = (int) numberOfEdges/20;
	fivePercentIncrements = fivePercentIncrements != 0 ? fivePercentIncrements : 1;

	for (int i=0; i<numberOfEdges; i++) {
		if (((*parameters).verbose) && ((i % fivePercentIncrements) == 0)) {
			int percentage = (i/fivePercentIncrements)*5;
			printf("%d%% done", percentage);
			if (percentage%20 == 0) {
				printf("\n");
			} else {
				printf(" •••• ");
			}
		}

		int fileFrom = 0, fileTo = 0;
		if (!fscanf(graphFile, "%d %d", &fileFrom, &fileTo)) {
			break;
		}

		apendElement(transitionMatrix, 1, fileFrom, fileTo);
	}

	// Calculates the outdegree of each page and assigns the uniform probability
	// of transition to the elements of the corresponding row
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
		double pageUniformProbability = 1. / pageOutdegree;
		for (int i=0; i<pageOutdegree; ++i) {
			if (currentElement->rowIndex == currentRow) {
				currentElement->value = pageUniformProbability;
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
	printf("%s [-c convergence_criterion] [-m max_iterations] [-a alpha] [-v] [-h] [-o output_filename] <graph_file>" \
		"\n-c convergence_criterion" \
		"\n\tthe convergence tolerance criterion" \
		"\n-m max_iterations" \
		"\n\tmaximum number of iterations to perform" \
		"\n-a alpha" \
		"\n\tthe damping factor" \
		"\n-v enable verbal output" \
		"\n-h enable history output to file" \
		"\n-o output_filename" \
		"\n\tfilename and path for the output" \
		"\n", programName);
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