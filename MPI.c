#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include "common/inc/timestamp.h"
#ifdef _OPENMP
#include <omp.h>
#else
int omp_get_thread_num(void) { return 0; }
int omp_get_num_threads(void) { return 1; }
#endif

int changed=0;// =1 exei ginei allagh , 0=dn exei ginei allagh
/* paral.c generations matrix_dimension */

void print(int ** temp,int mysize){
	int i,j;
	printf("\n--------------------------\n");
	for(i=0;i<=mysize+1;i++){
		if(i== mysize+1){printf("_________________________\n");}
		for(j=0;j<=mysize+1;j++){
			if(j==0 || j==mysize+1){
				printf(" |%d| ",temp[i][j]);
			}else{
				printf(" %d ",temp[i][j]);
			}
		}
		printf("\n");
		if(i==0 ){printf("_________________________\n");}
	}
	printf("\n--------------------------\n");
}

int neighbor_sum(int ** matrix,int i,int j){
	int sum=0;
	sum = matrix[i-1][j-1]+matrix[i-1][j]+
				matrix[i-1][j+1]+matrix[i][j-1]+matrix[i][j+1]+
				matrix[i+1][j-1]+matrix[i+1][j]+matrix[i+1][j+1];
	
	return sum;
}

int calculate_next_gen ( int number_of_neighbors ,int current_organism_state){
	if (number_of_neighbors < 2 || number_of_neighbors > 3)
		return 0;
	else if (number_of_neighbors == 3)
		return 1;
	else
		return current_organism_state;

}

double game_of_life(int matrix_size,int ntimes/*generations*/,MPI_Comm comm){
	
	
	int	rank, size ;/*rank= process id, size = number of processes*/
	int	right, left,up,down,up_right,up_left,down_right,down_left ;
	int	i, j, k;
	int	mysize, sum ;
	int	**matrix, **temp, **addr ;/* matrix einai to submatrix+2 diastash*/
	int	*left_send,*left_rcv,*right_send,*right_rcv; //auxiliary matrices gia na steiloume ta ka8eta stoixeia tou pinaka "orizontia" giAti den ta dexetai h send  ka8eta (h toul dn m rxetai twra tropos :P )
	double slavetime, totaltime, starttime ;
	int	my_offset;
	int sqrt_size; 
	
	/*Find my rank and the number of processes*/
	
	MPI_Comm_size(comm, &size) ;
	MPI_Comm_rank(comm, &rank) ;
	
	//Ypologismos ts tetragwnikhs rizas t plh8ous tn processes,k to size tou submatrix mas
	sqrt_size = sqrt(size);
	mysize = matrix_size/sqrt_size ;/*mysize= to size N tou upopinaka, matrix_size= to size t megalou pinaka,sqrt_size h tetragwnikh riza twn processes pou exoume*/
	
	
	/* Set neighbors */
	right = ((rank+1)%(sqrt_size))+sqrt_size*((int)rank/sqrt_size);
	left =  ((rank-1+sqrt_size)%(sqrt_size)+sqrt_size*((int)rank/sqrt_size)+size)%size;
	down =  (rank+sqrt_size)%size;
	up =  (rank+size-sqrt_size)%size;

	up_right = ((up+1)%(sqrt_size))+sqrt_size*((int)up/sqrt_size);
	up_left = ((up-1+sqrt_size)%(sqrt_size)+sqrt_size*((int)up/sqrt_size)+size)%size;
	down_right = ((down+1)%(sqrt_size))+sqrt_size*((int)down/sqrt_size);
	down_left = ((down-1+sqrt_size)%(sqrt_size)+sqrt_size*((int)down/sqrt_size)+size)%size;
	
	left_send = (int*)malloc(sizeof(int)*(mysize));
	left_rcv = (int*)malloc(sizeof(int)*(mysize));
	right_send = (int*)malloc(sizeof(int)*(mysize));
	right_rcv = (int*)malloc(sizeof(int)*(mysize));
	
	/*Allocate submatrix + temp gia to ouput + auxiliary matrices gia na steiloume
	 ta ka8eta stoixeia tou pinaka "orizontia"*/
	matrix = (int **)malloc(sizeof(int *)*(mysize+2)) ;
	temp = (int **)malloc(sizeof(int *)*(mysize+2)) ;
	for (i = 0; i < mysize+2; i++) {
		matrix[i] = (int *)malloc(sizeof(int)*(mysize+2));
		temp[i] = (int *)malloc(sizeof(int)*(mysize+2));
	}
	
	/* Initialize ta boundaries me 0  */
	for (j = 0; j < mysize+2; j++)
		matrix[0][j] = matrix[mysize+1][j] = temp[0][j] = temp[mysize+1][j] = 0;
	for (i = 0; i < mysize+2; i++)
		matrix[i][0] = matrix[i][mysize+1] = temp[i][0] = temp[i][mysize+1] = 0;
	
	/* Initialize to matrix randomly + initialize the auxiliary matrices */
	for (i = 1; i <= mysize; i++)  {
		for (j = 1; j<= mysize; j++){
			matrix[i][j] = (rand()%10==0)?1:0;
			temp[i][j] = 0;
			left_send[j] = 0;
			left_rcv[j] = 0;
			right_send[j] = 0;
			right_rcv[j] = 0;
		}
	}
	
	// xwnw ts ka8etes steiles tou submatrix m se temp orizonties gia na tis steilw
	for(i=1;i<mysize+1;i++){
		left_send[i-1] = matrix[i][1];
		right_send[i-1] = matrix[i][mysize];
	}
	
	
	
	/* 3ekiname to for gia tn upologismo twn generations*/
	starttime = MPI_Wtime();
	for(k = 0; k < ntimes; k++){
		/*Kanw initialize ta request k statuses */
		//printf("generation=%d",k);
		MPI_Request req[16];
		MPI_Status status[16];
	
	
		/* Send and receive boundary information */
		//Edw prp na kanw send receive
		//stelnw ston right,kanw receive apo ton right
		MPI_Isend(&right_send[0],mysize,MPI_INT,right,0,comm,req);//tou dinw thn de3ia sthlh
		MPI_Irecv(&right_rcv[0],mysize,MPI_INT,right,1,comm,req+1);//m dinei tn aristerh t ara tn de3ia m
		//stelnw ston left,kanw receive apo ton left
		MPI_Isend(&left_send[0],mysize,MPI_INT,left,1,comm,req+2);
		MPI_Irecv(&left_rcv[0],mysize,MPI_INT,left,0,comm,req+3);
		//stelnw ston up,kanw receive apo ton up
		MPI_Isend(&matrix[1][1],mysize,MPI_INT,up,2,comm,req+4);
		MPI_Irecv(&matrix[0][1],mysize,MPI_INT,up,3,comm,req+5);
		//stelnw ston down,kanw receive apo ton down
		MPI_Isend(&matrix[mysize][1],mysize,MPI_INT,down,3,comm,req+6);
		MPI_Irecv(&matrix[mysize+1][1],mysize,MPI_INT,down,2,comm,req+7);
		//stelnw ston up_left,kanw receive apo ton up_left
		MPI_Isend(&matrix[1][1],1,MPI_INT,up_left,5,comm,req+8);
		MPI_Irecv(&matrix[0][0],1,MPI_INT,up_left,6,comm,req+9);
		//stelnw ston down_right,kanw receive apo ton down_right
		MPI_Isend(&matrix[mysize][mysize],1,MPI_INT,down_right,6,comm,req+10);
		MPI_Irecv(&matrix[mysize+1][mysize+1],1,MPI_INT,down_right,5,comm,req+11);
		//stelnw ston down_left,kanw receive apo ton down_left
		MPI_Isend(&matrix[mysize][1],1,MPI_INT,down_left,7,comm,req+12);
		MPI_Irecv(&matrix[mysize+1][0],1,MPI_INT,down_left,8,comm,req+13);
		//stelnw ston up_right,kanw receive apo ton up_right
		MPI_Isend(&matrix[1][mysize],1,MPI_INT,up_right,8,comm,req+14);
		MPI_Irecv(&matrix[0][mysize+1],1,MPI_INT,up_right,7,comm,req+15);
	
		// KANW MIA FOR NA xwsw ta auxiliary left-right_rsv ston matrix m tn kanoniko
		for(i=0;i<mysize;i++){
			matrix[i+1][mysize+1] = right_rcv[i];
			matrix[i+1][0] = left_rcv[i];
		}
		//Ftiaxnw thn next generation gia to eswteriko p menei ametablhto apo to send rcv routine
		int n_sum=0;//neighbor_sum
		int has_changed=0;
//		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			int total = omp_get_num_threads();
			for(i=2; i<mysize ; i++){
//			#pragma omp for
				for(j=2; j<mysize; j++){
					n_sum = neighbor_sum(matrix,i,j);
					temp[i][j] = calculate_next_gen ( n_sum ,matrix[i][j]);
					if(matrix[i][j]!=temp[i][j]){//an exei diaforetikh timh o pinakas allazei eswterika dn paramenei idios
						has_changed=1;
					}
				}
			}
		}	

		//Kanw wait ola na staloun ta dedomena
		changed = 0;
		MPI_Waitall(16, req, status);
	
		//Ftiaxnw thn next generation gia to e3wteriko p ephreazotan apo to send rcv routine
		for (i = 1; i < mysize+1; i++){
			//upologizw to next gen gia tn de3ia akrh -sthlh tou eswterikou submatrix
			n_sum = neighbor_sum(matrix,i,mysize);
			temp[i][mysize] = calculate_next_gen (n_sum ,matrix[i][mysize]);
			//upologizw to next gen gia tn aristera akrh -sthlh tou eswterikou submatrix
			n_sum = neighbor_sum(matrix,i,1);
			temp[i][1] = calculate_next_gen (n_sum ,matrix[i][1]);
			//upologizw to next gen gia tn panw akrh -grammh tou eswterikou submatrix
			n_sum = neighbor_sum(matrix,1,i);
			temp[1][i] = calculate_next_gen (n_sum ,matrix[1][i]);
			//upologizw to next gen gia tn katw akrh -grammh tou eswterikou submatrix
			n_sum = neighbor_sum(matrix,mysize,i);
			temp[mysize][i] = calculate_next_gen (n_sum ,matrix[mysize][i]);
		}
		// Kanw swap ts pinakes
		addr = matrix;
		matrix = temp;
		temp = addr;
	
		
			/*Kanw all Reduce na dw an allazoun oi submatrices apo genia se genia*/
		MPI_Allreduce(&has_changed/*to sugkekrimeno process*/, &changed/*to sum*/, 1, MPI_INT, MPI_LOR,MPI_COMM_WORLD);
				//if(changed==0){
				//	printf("dn alla3an s auto to generation\n");
				//	break;
				//}
	
	}
	
	
}

// ./a.out generations dimension
int main (int argc,char *argv[]){
	
	int i,j;
	int generations,matrix_dimension;
	int rank, N, iters ;/*iters=generations*/
	double cpu;

	clock_t begin, end;
double time_spent;


/* here, do your time-consuming job */

//time_spent = (double)(end - begin) / CLOCKS_PER_SEC;


	//timestamp t_start;
	srand(1000);//(unsigned)time(NULL));
	
	if(argc<3){
		printf("lathos st argc\n");
		return 0;
	}
	iters = atoi(argv[1]);
	matrix_dimension = atoi(argv[2]);
	
	begin = clock();
	 
	/*Initialize MPI*/
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if(rank == 0)
		printf("gens %d dim %d\n", iters , matrix_dimension );
		
	/* Broadcast the size and # of iterations to all processes */
	MPI_Bcast(&matrix_dimension, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
	MPI_Bcast(&iters, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
	
	//t_start = getTimestamp();
	game_of_life( matrix_dimension, iters, MPI_COMM_WORLD );
	MPI_Finalize();

	end = clock();
	
	
	

	 
	
	
	
	return 0;
}
