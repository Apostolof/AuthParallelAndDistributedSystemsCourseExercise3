


//ALGORITHM 6 
//p_1 is period 1, p_2 is period 2 in no. of iterations 
float** filterMARP(int** A, float* x, float e, int max_iterations, float a, int p_1, int p_2, int size){
	
	float* x_new, x_c, y;
	int count = 0;
	float norm = 1.0;
	int* N, C, N_new, C_new;
	int** Ann, Acn;
	int i;
	/*********** ...malloc... *************/
	x_new = malloc(size*sizeof(float));
	
	x_c = malloc(size*sizeof(float));
	y = malloc(size*sizeof(float));
	N = malloc(size*sizeof(int));
	C = malloc(size*sizeof(int));
	N_new = malloc(size*sizeof(int));
	C_new = malloc(size*sizeof(int));
	Ann = malloc(size*sizeof(int *));
	for(i=0; i<size; i++){
		Ann[i] = malloc(size*sizeof(int));
	}
	Acn = malloc(size*sizeof(int *));
	for(i=0; i<size; i++){
		Acn[i] = malloc(size*sizeof(int));
	}
	
	/******************* start iterations ********************/
	while(count<max_iterations && norm<e){
		if(count == 0){
			gauss_Seidel(x, x_new, A, NULL, NULL, a, size);
		}
		else{
			gauss_Seidel(x, x_new, Ann, Acn, x_c, a, size);
		}
	
		count++;
		if(count%p_1==0){
			N = N_new; //might not be needed
			C = C_new;
			detect_Converged(N_new, C_new, x_new, x, e); //sinthiki stin prwti selida toy prwtou kefalaiou
			filterA(A, Ann, Acn, C_new, N_new, size);             //
			y =  compy(Acn, x, size);//Acn*x ;
			x_c = filterx(x_new, C, size);
		}
		if(count%p_2==0){
			norm = comp_norm(x_new, A); //||A*x - x||
		}
		x = x_new;
	}
	return x_new;
}

float* compy(int** Acn, float* x, int size){ //Acn*x  pollaplasiasmos
	float* y;
	y = malloc(size*sizeof(float));
	int i,j;
	for(i=0; i<size; i++){
		for(j=0; j<size; j++){
			y[i] += Acn[i][j]*x[j];
		}
	}
	return y;
}

void filterA(int** A, int** Ann, int** Acn, int* C_new, int* N_new, int size){
	int i,j;
	for(i=0; i<size; i++){
		if(C_new[i] == 1){
			for(j = 0; j<size;j++){
				Ann[i][j] = 0;    //sxesi 10 sto pdf
			}
		}
		else{
			for(j = 0; j<size;j++){
				Ann[i][j] = A[i][j];    //sxesi 10 sto pdf
			}
		}
	}
	
	
	//Acn part
}

float* filterx(float* xnew, int* C_new, int size){
	int i,j;
	float *x_c;
	x_c = malloc(size*sizeof(float));
	for(i=0; i<size; i++){
		if(C_new[i] == 1){
			x_c[i] = x[i];
		}
		else{
			x_c[i] = 0;
		}
	}
	return x_c;
}

void gauss_Seidel(float* x, float* x_new, int** Ann, int** Acn, float* x_c, float* y, float a, int size){
	int i, j;
	float sum1 = 0;
	float sum2 = 0;
	x_new = x; //giati ston algorithmo leei most updated values na xrisimopoiountai gia to sum1
	if(Acn == NULL && x_c == NULL){
		
		for(i=0; i<size; i++){
			for(j=0;j<i;j++){
				sum1 += Ann[i][j]*x_new[j]; //einai o pinakas A afou einai i prwti epanalipsi
			}
			for(j=i;j<size;j++){
				sum2 += Ann[i][j]*x[j];
			}
			x_new[i] = (1-a)+a*sum1+a*sum2;
		}
	}
	else{
		
		for(i=0; i<size; i++){
			for(j=0;j<i;j++){
				sum1 += Ann[i][j]*x_new[j];
			}
			for(j=i;j<size;j++){
				sum2 += Ann[i][j]*x[j];
			}
			
			x_new[i] = (1-a)+a*sum1+a*sum2+y[i]+x_c[i];
		}
		
	}

}

/*#include <stdio.h> 
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
		// break if end of file 
		if(m == EOF) break;
		// otherwise add one to the count of that particular character 
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

*/
