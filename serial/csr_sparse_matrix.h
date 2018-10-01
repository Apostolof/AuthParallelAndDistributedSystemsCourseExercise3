#ifndef CSR_SPARSE_MATRIX_H	/* Include guard */
#define CSR_SPARSE_MATRIX_H

/* ===== INCLUDES ===== */

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>

/* ===== STRUCTURES ===== */

// A sparse matrix in compressed SparseRow format.
typedef struct csrSparseMatrix {
	int size, numberOfNonZeroElements;
	int *rowCumulativeIndexes, *columnIndexes;
	double *values;
} CsrSparseMatrix;

/* ===== FUNCTION DEFINITIONS ===== */

// initCsrSparseMatrix creates and initializes the members of a CsrSparseMatrix
// structure instance.
CsrSparseMatrix initCsrSparseMatrix();

// allocMemoryForCsr allocates memory for the elements of the matrix.
void allocMemoryForCsr(CsrSparseMatrix *sparseMatrix, int numberOfElements);

// zeroOutRow assigns a zero value to all the elements of a row in the matrix.
void zeroOutRow(CsrSparseMatrix *sparseMatrix, int row);

// zeroOutColumn assigns a zero value to all the elements of a column in the
// matrix.
void zeroOutColumn(CsrSparseMatrix *sparseMatrix, int column);

// csrSparseMatrixVectorMultiplication calculates the product of a
// CsrSparseMatrix and a vector.
void csrSparseMatrixVectorMultiplication(CsrSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize);

// destroyCsrSparseMatrix frees all space used by the CsrSparseMatrix.
void destroyCsrSparseMatrix(CsrSparseMatrix *sparseMatrix);

// printCsrSparseMatrix prints the values of a CsrSparseMatrix.
void printCsrSparseMatrix(CsrSparseMatrix sparseMatrix);

#endif	// CSR_SPARSE_MATRIX_H