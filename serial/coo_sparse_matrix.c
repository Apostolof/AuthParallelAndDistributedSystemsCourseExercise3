#include "coo_sparse_matrix.h"

CooSparseMatrix initCooSparseMatrix() {
	CooSparseMatrix sparseMatrix;
	sparseMatrix.size = 0;
	sparseMatrix.elements = NULL;
	return sparseMatrix;
}

void allocMemoryForElements (CooSparseMatrix *sparseMatrix, int elements) {
	sparseMatrix->elements = (CooSparseMatrixElement **) malloc(
		elements * sizeof(CooSparseMatrixElement *));
}

void addElement(CooSparseMatrix *sparseMatrix, double value, int row, int column) {
	// Creates the new element
	CooSparseMatrixElement *newElement = (CooSparseMatrixElement *) malloc(
		sizeof(CooSparseMatrixElement));
	newElement->value = value;
	newElement->rowIndex = row;
	newElement->columnIndex = column;

	sparseMatrix->elements[sparseMatrix->size] = newElement;
	sparseMatrix->size = sparseMatrix->size + 1;
}

void zeroOutRow(CooSparseMatrix *sparseMatrix, int row) {
	for (int i=0; i<sparseMatrix->size; ++i) {
		CooSparseMatrixElement *element = sparseMatrix->elements[i];
		if (element->rowIndex == row) {
			element->value = 0;
		}
	}
}
void zeroOutColumn(CooSparseMatrix *sparseMatrix, int column) {
	for (int i=0; i<sparseMatrix->size; ++i) {
		CooSparseMatrixElement *element = sparseMatrix->elements[i];
		if (element->columnIndex == column) {
			element->value = 0;
		}
	}
}

int *getRowIndexes(CooSparseMatrix sparseMatrix, int row, int *rowSize) {
	*rowSize = 0;
	for (int i=0; i<sparseMatrix.size; ++i) {
		if (sparseMatrix.elements[i]->rowIndex == row) {
			++(*rowSize);
		}
	}

	if (!(*rowSize)) {
		return NULL;
	}

	int *indexes = (int *) malloc((*rowSize) * sizeof(int));
	int rowElementIndex = 0;
	for (int i=0; i<sparseMatrix.size; ++i) {
		if (sparseMatrix.elements[i]->rowIndex == row) {
			indexes[rowElementIndex] = i;
			++rowElementIndex;
		}
	}

	return indexes;
}

void transposeSparseMatrix(CooSparseMatrix *sparseMatrix) {
	for (int i=0; i<sparseMatrix->size; ++i) {
		CooSparseMatrixElement *element = sparseMatrix->elements[i];
		int tempRow = element->rowIndex;
		element->rowIndex = element->columnIndex;
		element->columnIndex = tempRow;
	}
}

void cooSparseMatrixVectorMultiplication(CooSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize) {
	// Initializes the elements of the product vector to zero
	for (int i=0; i<vectorSize; ++i) {
		(*product)[i] = 0;
	}

	CooSparseMatrixElement *element;
	for (int i=0; i<sparseMatrix.size; ++i) {
		element = sparseMatrix.elements[i];
		int row = element->rowIndex, column = element->columnIndex;

		if (row >= vectorSize) {
			printf("Error at sparseMatrixVectorMultiplication. Matrix has more rows than vector!\n");
			printf("row = %d\n", row);
			exit(EXIT_FAILURE);
		}

		(*product)[row] = (*product)[row] + element->value * vector[column];
	}
}

void destroyCooSparseMatrix(CooSparseMatrix *sparseMatrix) {
	for (int i=0; i<sparseMatrix->size; ++i) {
		free(sparseMatrix->elements[i]);
	}
	free(sparseMatrix->elements);
}

void printCooSparseMatrix(CooSparseMatrix sparseMatrix) {
	if (sparseMatrix.size == 0) {
		return;
	}

	CooSparseMatrixElement *element;
	for (int i=0; i<sparseMatrix.size; ++i) {
		element = sparseMatrix.elements[i];
		printf("[%d,%d] = %f\n", element->rowIndex, element->columnIndex,
			element->value);
	}
}