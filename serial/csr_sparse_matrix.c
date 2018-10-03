#include "csr_sparse_matrix.h"

CsrSparseMatrix initCsrSparseMatrix() {
	CsrSparseMatrix sparseMatrix;
	sparseMatrix.size = 0;
	sparseMatrix.numberOfElements = 0;

	sparseMatrix.values = NULL;
	sparseMatrix.columnIndexes = NULL;
	sparseMatrix.rowCumulativeIndexes = NULL;
	return sparseMatrix;
}

void allocMemoryForCsr(CsrSparseMatrix *sparseMatrix, int size, int numberOfElements) {
	sparseMatrix->values = (double *) malloc(numberOfElements * sizeof(double));
	sparseMatrix->columnIndexes = (int *) malloc(
		numberOfElements * sizeof(int));
	sparseMatrix->rowCumulativeIndexes = (int *) malloc(
		(size + 1) * sizeof(int));

	for (int i=0; i<size+1; ++i) {
		sparseMatrix->rowCumulativeIndexes[i] = 0;
	}

	sparseMatrix->size = size;
	sparseMatrix->numberOfElements = numberOfElements;
}

void zeroOutRow(CsrSparseMatrix *sparseMatrix, int row) {
	// Gets start and end indexes of the row's elements
	int startIndex = sparseMatrix->rowCumulativeIndexes[row],
	endIndex = sparseMatrix->rowCumulativeIndexes[row+1];
	for (int i=startIndex; i<endIndex; ++i) {
		sparseMatrix->values[i] = 0;
	}
}

void zeroOutColumn(CsrSparseMatrix *sparseMatrix, int column) {
	for (int i=0; i<sparseMatrix->numberOfElements; ++i){
		if(sparseMatrix->columnIndexes[i] == column){
			sparseMatrix->values[i] = 0;
		}
	}
}

void csrSparseMatrixVectorMultiplication(CsrSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize) {
	// Initializes the elements of the product vector to zero
	for (int i=0; i<vectorSize; ++i) {
		(*product)[i] = 0;
	}

	for (int i=0; i<sparseMatrix.size; ++i) {
		// Gets start and end indexes of this row's elements
		int startIndex = sparseMatrix.rowCumulativeIndexes[i],
		endIndex = sparseMatrix.rowCumulativeIndexes[i+1];

		if (startIndex == endIndex) {
			// This row has no elements
			continue;
		}

		double sum = 0;
		for(int j=startIndex; j<endIndex; ++j){
			int elementColumn = sparseMatrix.columnIndexes[j];
			sum += sparseMatrix.values[j] * vector[elementColumn];
		}

		(*product)[i] = sum;
	}
}

void destroyCsrSparseMatrix(CsrSparseMatrix *sparseMatrix) {
	free(sparseMatrix->values);
	free(sparseMatrix->rowCumulativeIndexes);
	free(sparseMatrix->columnIndexes);
}

void printCsrSparseMatrix(CsrSparseMatrix sparseMatrix) {
	if (sparseMatrix.size == 0) {
		return;
	}

	for (int i=0; i<sparseMatrix.size; ++i){
		int startIndex = sparseMatrix.rowCumulativeIndexes[i],
		endIndex = sparseMatrix.rowCumulativeIndexes[i+1];
		for(int j=startIndex; j<endIndex; ++j){
			printf("Row [%d] has [%d] nz elements: \n at column[%d] is value = %f \n",
				i, endIndex-startIndex,
				sparseMatrix.columnIndexes[j],
				sparseMatrix.values[j]);
		}
	}
}
