#include "sparse_matrix.h"

SparseMatrix createSparseMatrix() {
	SparseMatrix sparseMatrix;
	sparseMatrix.elements = 0;
	sparseMatrix.firstElement = NULL;
	sparseMatrix.lastElement = NULL;
	return sparseMatrix;
}

void apendElement(SparseMatrix *sparseMatrix, double value, int row, int column) {
	// Creates the new element
	SparseMatrixElement *newElement = (SparseMatrixElement *) malloc(sizeof(SparseMatrixElement));
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
		SparseMatrixElement *lastElement = sparseMatrix->lastElement;

		lastElement->nextElement = newElement;
		sparseMatrix->lastElement = newElement;
	}

	sparseMatrix->elements = sparseMatrix->elements + 1;
}

bool deleteElement(SparseMatrix *sparseMatrix, int row, int column) {
	if (sparseMatrix->elements == 0) {
		// Matrix is empty, nothing can be deleted
		return false;
	} else if (sparseMatrix->elements == 1) {
		// Matrix has one element. Deletes it.
		free(sparseMatrix->firstElement);
		sparseMatrix->firstElement = NULL;
		sparseMatrix->lastElement = NULL;
		sparseMatrix->elements = sparseMatrix->elements - 1;
		return true;
	}

	SparseMatrixElement *currentElement = sparseMatrix->firstElement;

	if (currentElement->rowIndex == row && currentElement->columnIndex == column) {
		sparseMatrix->firstElement = currentElement->nextElement;
		free(currentElement);
		sparseMatrix->elements = sparseMatrix->elements - 1;
		return true;
	}

	// Matrix has multiple elements. Finds the first element that has the coordinates
	// (row,column) and deletes it.
	for (int i=0; i<sparseMatrix->elements - 1; ++i) {
		SparseMatrixElement *nextElement = currentElement->nextElement;
		if (nextElement->rowIndex == row && nextElement->columnIndex == column) {
			currentElement->nextElement = nextElement->nextElement;
			if (currentElement->nextElement == NULL) {
				sparseMatrix->lastElement = currentElement;
			}
			free(nextElement);
			sparseMatrix->elements = sparseMatrix->elements - 1;
			return true;
		} else {
			currentElement = currentElement->nextElement;
		}
	}
}

SparseMatrixElement *getElement(SparseMatrix sparseMatrix, int row, int column) {
	SparseMatrixElement *currentElement = sparseMatrix.firstElement;
	do {
		if (currentElement->rowIndex == row && currentElement->columnIndex == column) {
			return currentElement;
		}
		currentElement = currentElement->nextElement;
	} while (currentElement != NULL);

	return NULL;
}

void transposeSparseMatrix(SparseMatrix *sparseMatrix) {
	SparseMatrixElement *currentElement = sparseMatrix->firstElement;
	for (int i=0; i<sparseMatrix->elements; ++i) {
		int temp = currentElement->rowIndex;
		currentElement->rowIndex = currentElement->columnIndex;
		currentElement->columnIndex = temp;

		currentElement = currentElement->nextElement;
	}
}

void sparseMatrixVectorMultiplication(SparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize) {
	// Initializes the elements of the product vector to zero
	for (int i=0; i<vectorSize; ++i) {
		(*product)[i] = 0;
	}

	SparseMatrixElement *element = sparseMatrix.firstElement;
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

void printSparseMatrix(SparseMatrix sparseMatrix) {
	if (sparseMatrix.elements == 0) {
		return;
	}

	SparseMatrixElement *currentElement = sparseMatrix.firstElement;
	for (int i=0; i<sparseMatrix.elements; ++i) {
		printf("[%d,%d] = %f\n", currentElement->rowIndex,
			currentElement->columnIndex, currentElement->value);
		currentElement = currentElement->nextElement;
	}
}