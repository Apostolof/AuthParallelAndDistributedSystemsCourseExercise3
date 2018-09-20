#ifndef COO_SPARSE_MATRIX_H	/* Include guard */
#define COO_SPARSE_MATRIX_H

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct cooSparseMatrixElement {
	double value;
	int rowIndex, columnIndex;
} CooSparseMatrixElement;

typedef struct cooSparseMatrix {
	int size;
	CooSparseMatrixElement **elements;
} CooSparseMatrix;

CooSparseMatrix initCooSparseMatrix();
void allocMemoryForElements (CooSparseMatrix *sparseMatrix, int elements);
void addElement(CooSparseMatrix *sparseMatrix, double value, int row, int column);
void zeroOutRow(CooSparseMatrix *sparseMatrix, int row);
void zeroOutColumn(CooSparseMatrix *sparseMatrix, int column);
int *getRowIndexes(CooSparseMatrix sparseMatrix, int row, int *rowSize);
void transposeSparseMatrix(CooSparseMatrix *sparseMatrix);
void cooSparseMatrixVectorMultiplication(CooSparseMatrix sparseMatrix, double *vector,
	double **product, int vectorSize);
void destroyCooSparseMatrix(CooSparseMatrix *sparseMatrix);
void printCooSparseMatrix(CooSparseMatrix sparseMatrix);

#endif	// COO_SPARSE_MATRIX_H