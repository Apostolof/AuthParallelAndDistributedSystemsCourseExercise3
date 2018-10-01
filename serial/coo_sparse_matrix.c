#include "coo_sparse_matrix.h"

CooSparseMatrix initCooSparseMatrix() {
	CooSparseMatrix sparseMatrix;
	sparseMatrix.size = 0;
	sparseMatrix.numberOfNonZeroElements = 0;
	sparseMatrix.elements = NULL;
	return sparseMatrix;
}

void allocMemoryForCoo(CooSparseMatrix *sparseMatrix, int numberOfElements) {
	sparseMatrix->elements = (CooSparseMatrixElement **) malloc(
		numberOfElements * sizeof(CooSparseMatrixElement *));
	sparseMatrix->size = numberOfElements;
}

void addElement(CooSparseMatrix *sparseMatrix, double value, int row, int column) {
	// Checks if there is enough space allocated
	if (sparseMatrix->numberOfNonZeroElements == sparseMatrix->size) {
		printf("Number of non zero elements exceeded size of matrix!\n");
		exit(EXIT_FAILURE);
	}

	// Creates the new element
	CooSparseMatrixElement *newElement = (CooSparseMatrixElement *) malloc(
		sizeof(CooSparseMatrixElement));
	newElement->value = value;
	newElement->rowIndex = row;
	newElement->columnIndex = column;

	// Adds the new element to the first empty (NULL) address of the matrix
	sparseMatrix->elements[sparseMatrix->numberOfNonZeroElements] = newElement;
	sparseMatrix->numberOfNonZeroElements = sparseMatrix->numberOfNonZeroElements + 1;
}

void transposeSparseMatrix(CooSparseMatrix *sparseMatrix) {
	for (int i=0; i<sparseMatrix->numberOfNonZeroElements; ++i) {
		CooSparseMatrixElement *element = sparseMatrix->elements[i];
		int tempRow = element->rowIndex;
		element->rowIndex = element->columnIndex;
		element->columnIndex = tempRow;
	}
}

/*
 * This function is a port of the one found here:
 * https://github.com/scipy/scipy/blob/3b36a57/scipy/sparse/sparsetools/coo.h#L34
*/
void transformToCSR(CooSparseMatrix initialSparseMatrix,
	CsrSparseMatrix *transformedSparseMatrix) {
	// Checks if the sizes of the two matrices fit
	if (initialSparseMatrix.numberOfNonZeroElements > transformedSparseMatrix->size) {
		printf("Transformed CSR matrix does not have enough space!\n");
		exit(EXIT_FAILURE);
	}

	// Calculates the elements per row
	for (int i=0; i<initialSparseMatrix.numberOfNonZeroElements; ++i){
		int rowIndex = initialSparseMatrix.elements[i]->rowIndex;
		transformedSparseMatrix->rowCumulativeIndexes[rowIndex] = 
		transformedSparseMatrix->rowCumulativeIndexes[rowIndex] + 1;
	}

	// Cumulative sums the non zero elements per row
	for (int i=0, sum=0; i<transformedSparseMatrix->size+1; ++i){
		int temp = transformedSparseMatrix->rowCumulativeIndexes[i];
		transformedSparseMatrix->rowCumulativeIndexes[i] = sum;
		sum += temp;
	}

	// Copies the values and columns of the elements
	for (int i=0; i<initialSparseMatrix.numberOfNonZeroElements; ++i){
		int row  = initialSparseMatrix.elements[i]->rowIndex;
		int destinationIndex = transformedSparseMatrix->rowCumulativeIndexes[row];

		transformedSparseMatrix->columnIndexes[destinationIndex] = initialSparseMatrix.elements[i]->columnIndex;
		transformedSparseMatrix->values[destinationIndex] = initialSparseMatrix.elements[i]->value;

		transformedSparseMatrix->rowCumulativeIndexes[row]++;
	}

	// Fixes the cumulative sum
	for (int i=0, last=0; i<=transformedSparseMatrix->size; i++){
		int temp = transformedSparseMatrix->rowCumulativeIndexes[i];
		transformedSparseMatrix->rowCumulativeIndexes[i] = last;
		last = temp;
	}

	transformedSparseMatrix->numberOfNonZeroElements = initialSparseMatrix.numberOfNonZeroElements;
}

void cooSparseMatrixVectorMultiplication(CooSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize) {
	// Initializes the elements of the product vector to zero
	for (int i=0; i<vectorSize; ++i) {
		(*product)[i] = 0;
	}

	CooSparseMatrixElement *element;
	for (int i=0; i<sparseMatrix.numberOfNonZeroElements; ++i) {
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
	for (int i=0; i<sparseMatrix->numberOfNonZeroElements; ++i) {
		free(sparseMatrix->elements[i]);
	}
	free(sparseMatrix->elements);
}

void printCooSparseMatrix(CooSparseMatrix sparseMatrix) {
	if (sparseMatrix.numberOfNonZeroElements == 0) {
		return;
	}

	CooSparseMatrixElement *element;
	for (int i=0; i<sparseMatrix.numberOfNonZeroElements; ++i) {
		element = sparseMatrix.elements[i];
		printf("[%d,%d] = %f\n", element->rowIndex, element->columnIndex,
			element->value);
	}
}