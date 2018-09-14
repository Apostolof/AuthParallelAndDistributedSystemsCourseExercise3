
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
			detect_Converged(N_new, C_new, x_new, x, e, size); //sinthiki stin prwti selida toy prwtou kefalaiou
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

void detect_Converged(int* N_new, int* C_new, float* x_new, float* x, float e, int size){
	int i;
	for( i=0; i<size; i++){
		if((norm(x_new[i] - x[i])/norm(x[i]))<e){ // 10 ^ -3
			//converged
			N_new[i] = 0;
			C_new[i] = 1;
		}
		else{
			N_new[i] = 1;
			C_new[i] = 0;
		}
	}
	return;
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
	
	
	//Acn part: pages that have converged to pages that haven't 
	// assuming A is column - based transition matrix (since it was transposed)
	
	for(i=0; i<size; i++){
		for(j=0; j<size ; j++){
			if(C_new[j] == 1 && N_new[i] == 1){
				Acn[i][j] = A[i][j];
			}
			else{
				Acn[i][j] = 0;
			}
		}
	}
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
