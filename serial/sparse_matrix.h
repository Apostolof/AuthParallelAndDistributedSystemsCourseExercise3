#ifndef SPARSE_MATRIX_H	/* Include guard */
#define SPARSE_MATRIX_H

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct sparseMatrixElement {
	double value;
	int rowIndex, columnIndex;
	struct sparseMatrixElement *nextElement;
} SparseMatrixElement;

typedef struct sparseMatrix {
	int elements;
	SparseMatrixElement *firstElement;
	SparseMatrixElement *lastElement;
} SparseMatrix;

SparseMatrix createSparseMatrix();
void apendElement(SparseMatrix *sparseMatrix, double value, int row, int column);
bool deleteElement(SparseMatrix *sparseMatrix, int row, int column);
SparseMatrixElement *getElement(SparseMatrix sparseMatrix, int row, int column);
void transposeSparseMatrix(SparseMatrix *sparseMatrix);
void sparseMatrixVectorMultiplication(SparseMatrix sparseMatrix, double *vector,
	double **product, int vectorSize);
void printSparseMatrix(SparseMatrix sparseMatrix);

#endif	// SPARSE_MATRIX_H