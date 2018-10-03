#include <sys/time.h>

#include "serial_gs_pagerank_functions.h"

struct timeval startwtime, endwtime;

int main(int argc, char **argv) {
	CsrSparseMatrix transitionMatrix = initCsrSparseMatrix();
	double *pagerankVector;
	bool convergenceStatus;
	Parameters parameters;
	int maxIterationsForConvergence;
	maxIterationsForConvergence = 0;
	parseArguments(argc, argv, &parameters);

	initialize(&transitionMatrix, &pagerankVector, &parameters);

	// Starts wall-clock timer
	gettimeofday (&startwtime, NULL);
	int* iterations = (int *)malloc(parameters.numberOfPages*sizeof(int));
	
	iterations = pagerank(&transitionMatrix, &pagerankVector,
		&convergenceStatus, parameters, &maxIterationsForConvergence);
	if (parameters.verbose) {
		printf(ANSI_COLOR_YELLOW "\n----- RESULTS -----\n" ANSI_COLOR_RESET);
		if (convergenceStatus) {
			printf(ANSI_COLOR_GREEN "Pagerank converged after %d iterations!\n" \
				ANSI_COLOR_RESET, maxIterationsForConvergence);
		} else {
			printf(ANSI_COLOR_RED "Pagerank did not converge after max number of" \
				" iterations (%d) was reached!\n" ANSI_COLOR_RESET, maxIterationsForConvergence);
		}
	}

	// Stops wall-clock timer
	gettimeofday (&endwtime, NULL);
	double seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 +
		endwtime.tv_sec - startwtime.tv_sec);
	printf("%s wall clock time = %f\n","Pagerank (Gauss-Seidel method), serial implementation",
		seq_time);
	if (!parameters.history) {
		// Always outputs numberOfPages, max_iterations, last pagerank and iterations 
		// for all pages
		savePagerankToFile(parameters.outputFilename, false, pagerankVector,
			parameters.numberOfPages, iterations, maxIterationsForConvergence);
	}

	free(pagerankVector);
	destroyCsrSparseMatrix(&transitionMatrix);
}
