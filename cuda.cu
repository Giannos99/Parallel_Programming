#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>
#include "./GoL.h"
#include "lcutil.h"
#include "timestamp.h"

__global__ void CycleRoutineGPU(char *CurrentState , char *NextState , int X , int Dimension){
	
	int sum=0;
	int id=blockIdx.x*blockDim.x + threadIdx.x;
	
	if (id<Dimension) {
		
		if(id+X<Dimension ) {
			sum += CurrentState[id+X];
		}
		if(id-X>=0){
			sum += CurrentState[id-X];
		}
		if(id/X == (id+1)/X) {
			sum += CurrentState[id+1];
		}
		if(id/X == (id-1)/X) {
			sum += CurrentState[id-1];
		}
		if(id+X<Dimension && (id+X)/X == (id+X+1)/X) {
			sum += CurrentState[id+X+1];
		}
		if(id+X<Dimension && (id+X)/X == (id+X-1)/X) {
			sum += CurrentState[id+X-1];
		}
		if(id-X>=0 && (id-X)/X == (id-X+1)/X) {
			sum += CurrentState[id-X+1];
		}
		if(id-X>=0 && (id-X)/X == (id-X-1)/X) {
			sum += CurrentState[id-X-1];
		}
		

		if (sum < 2 || sum > 3)
			NextState[id] = 0;
		else if (sum == 3)
			NextState[id] =  1;
		else
			NextState[id] = CurrentState[id];
					
	}
	
	__syncthreads();
}



void CycleRoutineCPU(char *CurrentState , char *NextState , int X , int Dimension){
	
	int sum=0;
	int id;
	for(id = 0 ; id < Dimension ; id++){
		if (id<Dimension) {
		
			if(id+X<Dimension ) {
				sum += CurrentState[id+X];
			}
			if(id-X>=0){
				sum += CurrentState[id-X];
			}
			if(id/X == (id+1)/X) {
				sum += CurrentState[id+1];
			}
			if(id/X == (id-1)/X) {
				sum += CurrentState[id-1];
			}
			if(id+X<Dimension && (id+X)/X == (id+X+1)/X) {
				sum += CurrentState[id+X+1];
			}
			if(id+X<Dimension && (id+X)/X == (id+X-1)/X) {
				sum += CurrentState[id+X-1];
			}
			if(id-X>=0 && (id-X)/X == (id-X+1)/X) {
				sum += CurrentState[id-X+1];
			}
			if(id-X>=0 && (id-X)/X == (id-X-1)/X) {
				sum += CurrentState[id-X-1];
			}
		

			if (sum < 2 || sum > 3)
				NextState[id] = 0;
			else if (sum == 3)
				NextState[id] =  1;
			else
				NextState[id] = CurrentState[id];
					
		}
	}
}


void CycleGPU(char *grid,int X , int Generations ){
	char *CurrentState,*NextState;
	int bytes=X*X*sizeof(char);
	int i;
	
	cudaMalloc((void**)&CurrentState,bytes);
	cudaMalloc((void**)&NextState,bytes);
	for(i=0;i<Generations;i++){
		cudaMemcpy( CurrentState , grid , bytes , cudaMemcpyHostToDevice );
		dim3 NumberOfThreads(X);
		dim3 NumberOfBlocks(X);
		CycleRoutineGPU<<<NumberOfBlocks,NumberOfThreads>>>(CurrentState,NextState,X,X*X);
		cudaMemcpy(grid,NextState,bytes,cudaMemcpyDeviceToHost);
	}
	cudaFree(NextState);
	cudaFree(CurrentState);
}

void CycleCPU(char *grid,int X , int Generations ){
	char *NextState;
	int bytes=X*X*sizeof(char);
	int i;
	
	//CurrentState = (char*)malloc(sizeof(char*)*X*X);
	NextState = (char*)malloc(sizeof(char*)*X*X);
	
	
	
	for(i=0;i<Generations;i++){
		CycleRoutineCPU(grid,NextState,X,X*X);
		memcpy( grid , NextState , bytes );
	}
}

int main(int argc , char **argv){
	int X,Generations,i,z;
	double cpu,gpu;
	timestamp t_start;
	X=atoi(argv[1]);
	Generations=atoi(argv[2]);
	
	char *grid_Gpu;
	char *grid_Cpu;

	grid_Gpu=(char*)malloc(sizeof(char*)*X*X);
	grid_Cpu=(char*)malloc(sizeof(char*)*X*X);
	
	srand(time(NULL));
	
	for(i=0;i<X;i++){
		for(z=0;z<X;z++){
			grid_Cpu[ i*X + z ] = grid_Gpu[ i*X + z ] = rand()%2;
		}
	}
	t_start = getTimestamp();
	CycleGPU(grid_Gpu , X, Generations );
	gpu = getElapsedtime(t_start);
	t_start = getTimestamp();
	CycleCPU(grid_Cpu , X, Generations );
	cpu = getElapsedtime(t_start);	
	free( grid_Gpu );
	free( grid_Cpu );
	printf("parallel: gens=%d dim=%d   ms %f   \n",Generations, X, gpu );
	printf("serial: gens=%d dim=%d   ms %f   \n\n",Generations, X, cpu );
	return 0;
}
