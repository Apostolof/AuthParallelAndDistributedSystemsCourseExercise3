#ifndef CSR_SPARSE_MATRIX_H	/* Include guard */
#define CSR_SPARSE_MATRIX_H

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct csrSparseMatrix {
	int size, numberOfNonZeroElements;
	int *rowCumulativeIndexes, *columnIndexes;
	double *values;
} CsrSparseMatrix;

CsrSparseMatrix initCsrSparseMatrix();
void allocMemoryForCsr(CsrSparseMatrix *sparseMatrix, int numberOfElements);
void zeroOutRow(CsrSparseMatrix *sparseMatrix, int row);
void zeroOutColumn(CsrSparseMatrix *sparseMatrix, int column);
void csrSparseMatrixVectorMultiplication(CsrSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize);
void destroyCsrSparseMatrix(CsrSparseMatrix *sparseMatrix);
void printCsrSparseMatrix(CsrSparseMatrix sparseMatrix);

#endif	// CSR_SPARSE_MATRIX_H