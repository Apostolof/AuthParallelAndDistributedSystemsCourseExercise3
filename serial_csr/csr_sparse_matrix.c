#include "csr_sparse_matrix.h"


CsrSparseMatrix initCsrSparseMatrix() {
	CsrSparseMatrix sparseMatrix;
	sparseMatrix.size = 0;   
	sparseMatrix.nnz = 0;  	 
	sparseMatrix.values = NULL;
	sparseMatrix.columnIndexes = NULL;
	sparseMatrix.rowaccInd = NULL;
	return sparseMatrix;
}

void allocMemoryForElements (CsrSparseMatrix *sparseMatrix, int edges) {
	sparseMatrix->values = (double *) malloc(
		edges * sizeof(double));				
	sparseMatrix->columnIndexes = (int *) malloc(
		edges * sizeof(int));
	sparseMatrix->rowaccInd = (int *) malloc(
		edges * sizeof(int));
}

void addElement(CsrSparseMatrix *sparseMatrix, double value, int row, int column) {
	
	for(int i=row+1; i<sparseMatrix->size; ++i){
		++sparseMatrix->rowaccInd[i];
	}
	sparseMatrix->nnz = sparseMatrix->nnz+1;
	sparseMatrix->values[sparseMatrix->rowaccInd[row-1]+1] = value;
	sparseMatrix->columnIndexes[sparseMatrix->rowaccInd[row-1]+1] = column;
	// Creates the new element
	/*CsrSparseMatrixElement *newElement = (CsrSparseMatrixElement *) malloc(
		sizeof(CsrSparseMatrixElement));
	
	newElement->value = value;
	newElement->rowIndex = row;
	newElement->columnIndex = column;

	sparseMatrix->elements[sparseMatrix->size] = newElement;   
	sparseMatrix->size = sparseMatrix->size + 1; */
}

void zeroOutRow(CsrSparseMatrix *sparseMatrix, int row) {
	int noofnnzinrow;
	if(row==0){
		noofnnzinrow = sparseMatrix->rowaccInd[row];
	}
	else{
		noofnnzinrow = sparseMatrix->rowaccInd[row]-sparseMatrix->rowaccInd[row-1];
	}
	int startdeleteInd = sparseMatrix->rowaccInd[row-1]+1;
	
	//delete the values and columnindexes of these rows by moving up the rest
	for(int i=0; i<noofnnzinrow; ++i){
		sparseMatrix->values[i+startdeleteInd] = sparseMatrix->values[sparseMatrix->nnz-noofnnzinrow+i];
		sparseMatrix->values[sparseMatrix->nnz-noofnnzinrow+i] = 0;
		sparseMatrix->columnIndexes[i+startdeleteInd] = sparseMatrix->columnIndexes[sparseMatrix->nnz-noofnnzinrow+i];
		sparseMatrix->columnIndexes[sparseMatrix->nnz-noofnnzinrow+i] = 0;
	}
	sparseMatrix->nnz = sparseMatrix->nnz - noofnnzinrow;
	
	//substract from accumulative no. of row nnz elements
	for(int i=row; i<sparseMatrix->size ; ++i){
		sparseMatrix->rowaccInd[i] -= noofnnzinrow;
	}
	
	/*for (int i=0; i<sparseMatrix->size; ++i) {
		CooSparseMatrixElement *element = sparseMatrix->elements[i];
		if (element->rowIndex == row) {
			element->value = 0;
		}
	}*/
}


void zeroOutColumn(CsrSparseMatrix *sparseMatrix, int column) {
	
	/*for (int i=0; i<sparseMatrix->size; ++i) {
		CooSparseMatrixElement *element = sparseMatrix->elements[i];
		if (element->columnIndex == column) {
			element->value = 0;
		}
	}
	*/
	for (int i=0; i<sparseMatrix->nnz; ++i){
		if(sparseMatrix->columnIndexes[i] == column){
			//delete columns by moving up the rest
			for(int j=i; j<sparseMatrix->nnz-1; ++j){
				sparseMatrix->columnIndexes[j] = sparseMatrix->columnIndexes[j+1];
				sparseMatrix->values[j] = sparseMatrix->values[j+1];
			}
			int flag = 0;
			//adjust rowaccInd 
			for(int j=0; j<sparseMatrix->size; ++j){
				if(sparseMatrix->rowaccInd[j] > i){
					flag = 1;        //must be substracted since column belonged to this row
				}
				if(flag){
					--sparseMatrix->rowaccInd[j];  //substract till end of rows 
				}
			}
			
		}
	}
	
}

int *getRowIndexes(CsrSparseMatrix sparseMatrix, int row, int *rowSize) { 
	*rowSize = 0;
	/*for (int i=0; i<sparseMatrix.size; ++i) {
		if (sparseMatrix.elements[i]->rowIndex == row) {
			++(*rowSize);
		}
	}

	if (!(*rowSize)) {
		return NULL;
	}*/
	if((row-1)>0 && (sparseMatrix.rowaccInd[row]-sparseMatrix.rowaccInd[row-1])>0){
		(*rowSize) = sparseMatrix.rowaccInd[row]-sparseMatrix.rowaccInd[row-1];
	}
	else if((sparseMatrix.rowaccInd[row]-sparseMatrix.rowaccInd[row-1])>0){ //if row = 0
		(*rowSize) = sparseMatrix.rowaccInd[row];
	}
	else{
		return NULL;
	}
 
	int *indexes = (int *) malloc((*rowSize) * sizeof(int));
	for (int i=1; i<=(*rowSize); ++i) {
		
		indexes[i-1] = sparseMatrix.rowaccInd[row-1]+i;				
		
	}

	return indexes;
}

void transposeSparseMatrix(CsrSparseMatrix *sparseMatrix) {
	/*for (int i=0; i<sparseMatrix->size; ++i) {
		CooSparseMatrixElement *element = sparseMatrix->elements[i];
		int tempRow = element->rowIndex;
		element->rowIndex = element->columnIndex;
		element->columnIndex = tempRow;
	}*/
	double* values_t = (double *) malloc(
		sparseMatrix->size * sizeof(double));
	int* rowIndexes = (int *) malloc(
		sparseMatrix->size * sizeof(int));
	int* colaccInd = (int *) malloc(
		sparseMatrix->size * sizeof(int));
		
	
	
	int columncount, nnznew = 0;
	//for all columns
	for(columncount = 0; columncount<sparseMatrix->size; ++columncount){
		//index for searching in columnIndexes matrix
		for(int i = 0; i<sparseMatrix->nnz;++i){
			if(sparseMatrix->columnIndexes[i] == columncount){
				//Find which row it belongs to
				for(int j=0; j<sparseMatrix->size; ++j){
					if(sparseMatrix->rowaccInd[j] == i){
						rowIndexes[nnznew] = j-1;
						values_t[nnznew] = sparseMatrix->values[i];
						for(int k=i; k<sparseMatrix->size; ++k){
							++colaccInd[k];
						}
						++nnznew;
					}
				}
				
			}
		}
	}
	
	memcpy(sparseMatrix->values, values_t, sparseMatrix->size*sizeof(double));
	memcpy(sparseMatrix->columnIndexes, rowIndexes, sparseMatrix->size*sizeof(int));
	memcpy(sparseMatrix->rowaccInd, colaccInd, sparseMatrix->size*sizeof(int) );
	sparseMatrix->nnz = nnznew;
}

void csrSparseMatrixVectorMultiplication(CsrSparseMatrix sparseMatrix,
	double *vector, double **product, int vectorSize) {
	// Initializes the elements of the product vector to zero
	for (int i=0; i<vectorSize; ++i) {
		(*product)[i] = 0;
	}

	/*CooSparseMatrixElement *element;
	for (int i=0; i<sparseMatrix.size; ++i) {
		element = sparseMatrix.elements[i];
		int row = element->rowIndex, column = element->columnIndex;

		if (row >= vectorSize) {
			printf("Error at sparseMatrixVectorMultiplication. Matrix has more rows than vector!\n");
			printf("row = %d\n", row);
			exit(EXIT_FAILURE);
		}

		(*product)[row] = (*product)[row] + element->value * vector[column];
	}*/
	int t;
	//for every row
	for (int i=0; i<sparseMatrix.size; ++i) {
		if(i==0){
			t = sparseMatrix.rowaccInd[0];
		}
		else{
			t = sparseMatrix.rowaccInd[i]-sparseMatrix.rowaccInd[i-1];
		}
		for(int j=0; j<t; ++j){
			for(int k=0; k<vectorSize; ++k){
				if(sparseMatrix.columnIndexes[sparseMatrix.rowaccInd[i]+t]==k){
					(*product)[k] += sparseMatrix.values[sparseMatrix.rowaccInd[i]+t]*vector[k];
				}
				else if(sparseMatrix.columnIndexes[sparseMatrix.rowaccInd[i]+t]>k){
					printf("Error at sparseMatrixVectorMultiplication. Matrix has more columns than vector rows!\n");
					exit(EXIT_FAILURE);
				}
			}
			
		}
	}
	
	
}

void destroyCsrSparseMatrix(CsrSparseMatrix *sparseMatrix) {
	/*for (int i=0; i<sparseMatrix->size; ++i) {
		free(sparseMatrix->elements[i]);
	}*/
	free(sparseMatrix->values);
	free(sparseMatrix->rowaccInd);
	free(sparseMatrix->columnIndexes);
 
}

void printCsrSparseMatrix(CsrSparseMatrix sparseMatrix) {
	if (sparseMatrix.size == 0) {
		return;
	}
	/*
	CooSparseMatrixElement *element;
	for (int i=0; i<sparseMatrix.size; ++i) {
		element = sparseMatrix.elements[i];
		printf("[%d,%d] = %f\n", element->rowIndex, element->columnIndex,
			element->value);
	}*/
	int t;
	for (int i=0; i<sparseMatrix.size; ++i){
		if(i==0){
			t = sparseMatrix.rowaccInd[i];
		}
		else{
			t = sparseMatrix.rowaccInd[i]-sparseMatrix.rowaccInd[i-1];
		}
		for(int j=0; j<t ; ++j){
			printf("Row [%d] has [%d] nz elements: \n at column[%d] is value = %f \n",
				i, t, sparseMatrix.columnIndexes[sparseMatrix.rowaccInd[i]+j], sparseMatrix.values[sparseMatrix.rowaccInd[i]+j]);
		}
	}
}