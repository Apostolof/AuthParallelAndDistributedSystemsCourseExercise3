#include <sys/time.h>

#include "serial_gs_pagerank_functions.h"

struct timeval startwtime, endwtime;
double seq_time;

int main(int argc, char **argv) {
	int **directedWebGraph;
	double **transitionMatrix, *pagerankVector;
	Parameters parameters;

	parseArguments(argc, argv, &parameters);

	initialize(&directedWebGraph, &transitionMatrix, &pagerankVector, &parameters);

	// Starts wall-clock timer
	gettimeofday (&startwtime, NULL);

	int iterations = pagerank(&transitionMatrix, &pagerankVector, parameters);
	if (parameters.verbose) {
		printf("\n----- Results -----\
			\nTotal iterations = %d\n", iterations);
	}

	// Stops wall-clock timer
	gettimeofday (&endwtime, NULL);
	double seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 +
		endwtime.tv_sec - startwtime.tv_sec);
	printf("%s wall clock time = %f\n","Pagerank (Gauss-Seidel method), serial implementation",
		seq_time);
}
