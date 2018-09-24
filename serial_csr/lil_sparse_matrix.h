#ifndef LIL_SPARSE_MATRIX_H	/* Include guard */
#define LIL_SPARSE_MATRIX_H

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct lilSparseMatrixElement {
	double value;
	int rowIndex, columnIndex;
	struct lilSparseMatrixElement *nextElement;
} LilSparseMatrixElement;

typedef struct lilSparseMatrix {
	int elements;
	LilSparseMatrixElement *firstElement;
	LilSparseMatrixElement *lastElement;
} LilSparseMatrix;

LilSparseMatrix createLilSparseMatrix();
void apendElement(LilSparseMatrix *sparseMatrix, double value, int row,
	int column);
void lilSparseMatrixVectorMultiplication(LilSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize);
void destroyLilSparseMatrix(LilSparseMatrix *sparseMatrix);
void printLilSparseMatrix(LilSparseMatrix sparseMatrix);

#endif	// LIL_SPARSE_MATRIX_H