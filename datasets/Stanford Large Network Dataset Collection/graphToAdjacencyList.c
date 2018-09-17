#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct pair {
	int from;
	int to;
} Pair;

int main(int argc, char const *argv[]) {
	if (argc != 2) {
		printf("Usage: ./graphToAdjacencyList <graphFile>\n");
		exit(EXIT_FAILURE);
	}

	FILE *graphFile = fopen(argv[1], "r");
	if (!graphFile) {
		printf("Error opening file \n");
		exit(EXIT_FAILURE);
	}

	char buffer[512];
	// Skips the first two lines
	fgets(buffer, 512, graphFile);
	fgets(buffer, 512, graphFile);

	// Third line has the numbers of nodes and edges, we need to parse them.
	int numberOfNodes = 0, numberOfEdges;
	fgets(buffer, 512, graphFile);

	char *token = strtok(buffer, " ");
	int getNodes = 0, getEdges = 0;

	while (token != NULL) {
		if (strcmp(token, "Nodes:") == 0) {
			getNodes = 1;
		} else if (getNodes == 1) {
			numberOfNodes = atoi(token);
			getNodes = 0;
		} else if (strcmp(token, "Edges:") == 0) {
			getEdges = 1;
		} else if (getEdges == 1) {
			numberOfEdges = atoi(token);
			break;
		}
		token = strtok (NULL, " ,.-");
	}

	// Skips the fourth line
	fgets(buffer, 512, graphFile);

	Pair **edges = (Pair **) malloc(numberOfEdges * sizeof(Pair *));
	printf("Reading edges from file...\n");
	for (int i=0; i<numberOfEdges; i++) {
		edges[i] = (Pair *) malloc(sizeof(Pair));

		int f_from = 0, f_to = 0;
		if (!fscanf(graphFile, "%d %d", &f_from, &f_to)) {
			break;
		}

		edges[i]->from = f_from;
		edges[i]->to = f_to;
	}


	FILE *adjacentListFile = fopen("adjacentListCreated", "w");
	int index = 0;

	printf("\nWriting nodes to file...\n");
	for (int i=0; i<numberOfNodes; ++i) {
		int hasOutlinks = 0;

		fprintf(adjacentListFile, "%d: ", i);

		while (index < numberOfEdges && edges[index]->from == i) {
			fprintf(adjacentListFile, "%d ", edges[index]->to);
			if (!hasOutlinks) hasOutlinks = 1;
			++index;
		}

		if (!hasOutlinks) {
			fprintf(adjacentListFile, "-1 ");
		}

		fprintf(adjacentListFile, "\n");
	}

	return 0;
}