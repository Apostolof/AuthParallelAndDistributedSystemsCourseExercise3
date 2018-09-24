#ifndef CSR_SPARSE_MATRIX_H	/* Include guard */
#define CSR_SPARSE_MATRIX_H

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct csrSparseMatrix {
	double* values;
	int* rowaccInd;   				//without the first cell, always 0
	int* columnIndexes;
	int size;						//no. of rows
	int nnz;						//no. of non zero elements

} CsrSparseMatrix;

CsrSparseMatrix initCsrSparseMatrix();
void allocMemoryForElements (CsrSparseMatrix *sparseMatrix, int edges);
void addElement(CsrSparseMatrix *sparseMatrix, double value, int row, int column);
void zeroOutRow(CsrSparseMatrix *sparseMatrix, int row);
void zeroOutColumn(CsrSparseMatrix *sparseMatrix, int column);
int *getRowIndexes(CsrSparseMatrix sparseMatrix, int row, int *rowSize);
void transposeSparseMatrix(CsrSparseMatrix *sparseMatrix);
void csrSparseMatrixVectorMultiplication(CsrSparseMatrix sparseMatrix, double *vector,
	double **product, int vectorSize);
void destroyCsrSparseMatrix(CsrSparseMatrix *sparseMatrix);
void printCsrSparseMatrix(CsrSparseMatrix sparseMatrix);

#endif	// CSR_SPARSE_MATRIX_H