#include "lil_sparse_matrix.h"

LilSparseMatrix createLilSparseMatrix() {
	LilSparseMatrix sparseMatrix;
	sparseMatrix.elements = 0;
	sparseMatrix.firstElement = NULL;
	sparseMatrix.lastElement = NULL;
	return sparseMatrix;
}

void apendElement(LilSparseMatrix *sparseMatrix, double value, int row, int column) {
	// Creates the new element
	LilSparseMatrixElement *newElement = (LilSparseMatrixElement *) malloc(sizeof(LilSparseMatrixElement));
	newElement->value = value;
	newElement->rowIndex = row;
	newElement->columnIndex = column;
	newElement->nextElement = NULL;

	if (sparseMatrix->firstElement == NULL) {
		// Sparse matrix is empty, this is the first element
		sparseMatrix->firstElement = newElement;
		sparseMatrix->lastElement = newElement;
	} else {
		//Gets last element of the matrix
		LilSparseMatrixElement *lastElement = sparseMatrix->lastElement;

		lastElement->nextElement = newElement;
		sparseMatrix->lastElement = newElement;
	}

	sparseMatrix->elements = sparseMatrix->elements + 1;
}

void lilSparseMatrixVectorMultiplication(LilSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize) {
	// Initializes the elements of the product vector to zero
	for (int i=0; i<vectorSize; ++i) {
		(*product)[i] = 0;
	}

	LilSparseMatrixElement *element = sparseMatrix.firstElement;
	for (int i=0; i<sparseMatrix.elements; ++i) {
		int row = element->rowIndex, column = element->columnIndex;

		if (row >= vectorSize) {
			printf("Error at sparseMatrixVectorMultiplication. Matrix has more rows than vector!\n");
			printf("row = %d\n", row);
			exit(EXIT_FAILURE);
		}

		(*product)[row] = (*product)[row] + element->value * vector[column];
		element = element->nextElement;
	}
}

void destroyLilSparseMatrix(LilSparseMatrix *sparseMatrix) {
	LilSparseMatrixElement *currentElement = sparseMatrix->firstElement;
	while (currentElement != NULL) {
		LilSparseMatrixElement *toDelete = currentElement;
		currentElement = currentElement->nextElement;
		free(toDelete);
	}
}

void printLilSparseMatrix(LilSparseMatrix sparseMatrix) {
	if (sparseMatrix.elements == 0) {
		return;
	}

	LilSparseMatrixElement *currentElement = sparseMatrix.firstElement;
	for (int i=0; i<sparseMatrix.elements; ++i) {
		printf("[%d,%d] = %f\n", currentElement->rowIndex,
			currentElement->columnIndex, currentElement->value);
		currentElement = currentElement->nextElement;
	}
}