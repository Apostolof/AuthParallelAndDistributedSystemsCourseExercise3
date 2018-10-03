#include <sys/time.h>

#include "openmp_gs_pagerank_functions.h"

struct timeval startwtime, endwtime;

int main(int argc, char **argv) {
	CsrSparseMatrix transitionMatrix = initCsrSparseMatrix();
	double *pagerankVector;
	bool convergenceStatus;
	Parameters parameters;
	int maxIterationsForConvergence = 0;

	parseArguments(argc, argv, &parameters);

	initialize(&transitionMatrix, &pagerankVector, &parameters);

	// Saves information about the dataset to the output file
	{
		FILE *outputFile;
		outputFile = fopen(parameters.outputFilename, "w");

		if (outputFile == NULL) {
			printf("Error while opening the output file.\n");
			exit(EXIT_FAILURE);
		}

		fprintf(outputFile, "Pagerank will run for the dataset %s\n"\
			"Dataset contains %d pages with %d outlinks.\n",
			parameters.graphFilename, parameters.numberOfPages, transitionMatrix.size);

		fclose(outputFile);
	}

	// Starts wall-clock timer
	gettimeofday (&startwtime, NULL);
	int* iterations = (int *)malloc(parameters.numberOfPages*sizeof(int));

	// Calculates pagerank
	iterations = pagerank(&transitionMatrix, &pagerankVector,
		&convergenceStatus, parameters, &maxIterationsForConvergence);

	// Stops wall-clock timer
	gettimeofday (&endwtime, NULL);
	double seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 +
		endwtime.tv_sec - startwtime.tv_sec);
	printf("%s wall clock time = %f\n","Pagerank (Gauss-Seidel method), serial implementation",
		seq_time);

	printf(ANSI_COLOR_YELLOW "\n----- RESULTS -----\n" ANSI_COLOR_RESET);
	if (convergenceStatus) {
		printf(ANSI_COLOR_GREEN "Pagerank converged after %d iterations!\n" \
			ANSI_COLOR_RESET, maxIterationsForConvergence);
	} else {
		printf(ANSI_COLOR_RED "Pagerank did not converge after max number of" \
			" iterations (%d) was reached!\n" ANSI_COLOR_RESET, maxIterationsForConvergence);
	}

	// Saves results to the output file
	savePagerankToFile(parameters.outputFilename, iterations, pagerankVector,
		parameters.numberOfPages, maxIterationsForConvergence);

	free(pagerankVector);
	destroyCsrSparseMatrix(&transitionMatrix);
}
