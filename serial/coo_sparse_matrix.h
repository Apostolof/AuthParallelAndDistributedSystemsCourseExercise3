#ifndef COO_SPARSE_MATRIX_H	/* Include guard */
#define COO_SPARSE_MATRIX_H

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>

#include "csr_sparse_matrix.h"

typedef struct cooSparseMatrixElement {
	double value;
	int rowIndex, columnIndex;
} CooSparseMatrixElement;

typedef struct cooSparseMatrix {
	int size, numberOfNonZeroElements;
	CooSparseMatrixElement **elements;
} CooSparseMatrix;

CooSparseMatrix initCooSparseMatrix();
void allocMemoryForCoo(CooSparseMatrix *sparseMatrix, int numberOfElements);
void addElement(CooSparseMatrix *sparseMatrix, double value, int row,
	int column);
void transposeSparseMatrix(CooSparseMatrix *sparseMatrix);
void transformToCSR(CooSparseMatrix initialSparseMatrix,
	CsrSparseMatrix *transformedSparseMatrix);
void cooSparseMatrixVectorMultiplication(CooSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize);
void destroyCooSparseMatrix(CooSparseMatrix *sparseMatrix);
void printCooSparseMatrix(CooSparseMatrix sparseMatrix);

#endif	// COO_SPARSE_MATRIX_H