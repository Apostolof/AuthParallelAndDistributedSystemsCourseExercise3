#ifndef COO_SPARSE_MATRIX_H	/* Include guard */
#define COO_SPARSE_MATRIX_H

/* ===== INCLUDES ===== */

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>

#include "csr_sparse_matrix.h"

/* ===== STRUCTURES ===== */

// One element of the coordinate formated sparse matrix.
typedef struct cooSparseMatrixElement {
	double value;
	int rowIndex, columnIndex;
} CooSparseMatrixElement;

// A sparse matrix in COOrdinate format (aka triplet format).
typedef struct cooSparseMatrix {
	int size, numberOfNonZeroElements;
	CooSparseMatrixElement **elements;
} CooSparseMatrix;

/* ===== FUNCTION DEFINITIONS ===== */

// initCooSparseMatrix creates and initializes the members of a CooSparseMatrix
// structure instance.
CooSparseMatrix initCooSparseMatrix();

//allocMemoryForCoo allocates memory for the elements of the matrix.
void allocMemoryForCoo(CooSparseMatrix *sparseMatrix, int numberOfElements);

// addElement adds an element representing the triplet passed in the arguments
// to the first empty address of the space allocated for the elements.
void addElement(CooSparseMatrix *sparseMatrix, double value, int row,
	int column);

// transposeSparseMatrix transposes the matrix.
void transposeSparseMatrix(CooSparseMatrix *sparseMatrix);

// transformToCSR transforms the sparse matrix representation format from COO
// to CSR.
void transformToCSR(CooSparseMatrix initialSparseMatrix,
	CsrSparseMatrix *transformedSparseMatrix);

// cooSparseMatrixVectorMultiplication calculates the product of a
// CooSparseMatrix and a vector.
void cooSparseMatrixVectorMultiplication(CooSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize);

// destroyCooSparseMatrix frees all space used by the CooSparseMatrix.
void destroyCooSparseMatrix(CooSparseMatrix *sparseMatrix);

// printCooSparseMatrix prints the values of a CooSparseMatrix.
void printCooSparseMatrix(CooSparseMatrix sparseMatrix);

#endif	// COO_SPARSE_MATRIX_H