#include <stdio.h> 
#include <stdlib.h>


void main(){
	printf("Hello I am maria \n");
	
	//Read from file of adjacency matrix
	FILE *adjm;
	
	int x;
	adjm = fopen("adj_matrix.txt", "r+");
	if(!adjm){
		printf("Error opening file \n");
		exit(0);
	}
	
	//Read dimensions of file (square)
	int m, count=0;
	while((m=fgetc(adjm))) {
		/* break if end of file */
		if(m == EOF) break;
		/* otherwise add one to the count of that particular character */
		if(m=='\n'){
			count+=1;
		}
    }
	printf("Line count of matrix is %d \n", count);
	int d = count;
	int i,j;
	
	// Put values in matrix A
	int** A = malloc(d*sizeof(int *));
	for( i=0; i<d ; i++){
		A[i] = malloc(d*sizeof(int));
	}
	for( i=0; i<d ; i++){
		for(j=0 ; j<d; j++){
			if(!fscanf(adjm, "%d ", &A[i][j])){
				break;
			}
			//printf("A[%d][%d] = %d", i , j, A[i][j]);
		}
	}
	fclose(adjm);
	printf(" First val is %d Last val is %d \n", A[0][0], A[d-1][d-1]);
	
	//Make A appropriate for the algorithm
	// no page has outdegree 0, using uniform probability 1/n, no personalization
	
	int* flag;
	flag = malloc(d*sizeof(int));
	for(i=0; i<d ; i++){
		flag[i] = 0;
	}
	for( i=0; i<d ; i++){
		for(j=0 ; j<d; j++){
			if(A[i][j]!=0){
				flag[i]=1;
			}
		}
		if(flag[i] == 1){
			for(j=0 ; j<d; j++){
				A[i][j] = 1;     
			}
		}
		printf("A[%d][%d] = %d", i , j, A[i][j]);
	}
	
	//Change to transpose of matrix
	// Rows become columns
	
	int **AT = malloc(d*sizeof(int *));
	for( i=0; i<d ; i++){
		AT[i] = malloc(d*sizeof(int));
	}
	
	for( i=0; i<d ; i++){
		for(j=0 ; j<d; j++){
			AT[j][i] = A[i][j];
		}
	}
	
	



}