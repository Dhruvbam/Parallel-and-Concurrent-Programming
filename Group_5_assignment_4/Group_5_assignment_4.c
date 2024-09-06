#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/*
Group 5
Brandon Bratcher    R#11520450
James Gibson        R#11671129
Richard Pearson     R#11786900
*/


void default_arr(int, int, int**);
void hw4(int, int**, int**);
void ij_loop1(int,int**,int**,int*,int*);
void ij_loop2(int,int**,int**);
void get_k_col_and_row(int, int, int**, int*, int*);
int rowwise_one_to_many(int*,int,int,int);
int colwise_one_to_many(int*,int,int,int);
int isSquare(int);
int find_min(int,int,int);

int main(int argc, char **argv) {
    //variable declaration
    int pid, np;
    //MPI variable initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //gaurantee only a square number of processes are provided
    if (!isSquare(np)) {
        return -1;
    }

    //size variables
    int sqrt_np = sqrt(np);
    int chunk_size = 3;
    int n = sqrt_np*chunk_size;
    
    //array 1
    int **arr = (int**)malloc(chunk_size * sizeof(int*));
    for (int i=0; i<chunk_size; i++) {
        arr[i] = (int*)malloc(chunk_size * sizeof(int));
    }

    //array 1 generation
    default_arr(n, sqrt_np, arr);
    //print_2d_ptr(n, n, arr);
    
    //array 2
    int **arr2 = (int**)malloc(n * sizeof(int*));
    for (int i=0; i<n; i++) {
        arr2[i] = (int*)malloc(n * sizeof(int));
    }

    //display array 1
    printf("D0 pid:%d\n", pid);
    for (int i=0; i<chunk_size; i++) {
        for (int j=0; j<chunk_size; j++) {
            printf("%d ", arr[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    //perform algorithm on array 1 to generate arrat 2
    hw4(n, arr, arr2);

    //display array 2
    printf("D pid:%d\n", pid);
    for (int i=0; i<chunk_size; i++) {
        for (int j=0; j<chunk_size; j++) {
            printf("%d ", arr2[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    
    
    //exit
    MPI_Finalize();
    return 0;
}

//array generation function (allocates each section to appropriate process)
void default_arr(int n, int sqrt_np, int** arr) {
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    int chunk_size = n/sqrt_np;
    
    for (int i = 0; i<n; i++) {
        int row = i/chunk_size;    
        for (int j=0; j<n; j++) {
            int col = j/chunk_size;
            int curr_process = row*sqrt_np + col;
            
            if (pid == curr_process) {
                arr[i%chunk_size][j%chunk_size] = 1;
            }
        }
        int curr_process = row*sqrt_np + row;
        if (pid == curr_process) {
            arr[i%chunk_size][i%chunk_size] = 0;
        }
    }
}

//d0->D algorithm
void hw4(int n, int** D0, int** D) {
    //mpi variables
    int pid, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //size variables
    int sqrt_np = sqrt(np);
    int chunk_size = n/sqrt_np;
    
    //arrays for dynamically updated and loaded data
    int *k_col = (int*)malloc(sqrt_np * sizeof(int));
    int *k_row = (int*)malloc(sqrt_np * sizeof(int));
    
    //processor id's holding the dynamically updated value needed
    int row_sender = (pid/sqrt_np)*sqrt_np;
    int col_sender = pid%sqrt_np;
    
    //k loop
    for (int k=0; k<n; k++) {
        //recalculate processor id's when necessary
        if (k%chunk_size == 0 && k!=0) {
            col_sender += 1;
            row_sender += sqrt_np;
        }
        
        //send or receive dynamically updated arrays
        get_k_col_and_row(n,k,D0,k_col,k_row);
        rowwise_one_to_many(k_col, chunk_size, k, row_sender);
        colwise_one_to_many(k_row, chunk_size, k, col_sender);

        //run parallelized loops
        ij_loop1(chunk_size,D0,D,k_col,k_row);
        ij_loop2(chunk_size,D0,D);
    }
}

//parallelized loop 1
void ij_loop1(int chunk_size, int** D0, int** D, int* k_col, int* k_row) {
    for (int i=0; i<chunk_size; i++) {
        for (int j=0; j<chunk_size; j++) {
            D[i][j] = find_min(D0[i][j], k_col[i], k_row[j]);
        }
    }
}

//parallelized loop 2
void ij_loop2(int chunk_size, int** D0, int** D) {
    for (int i=0; i<chunk_size; i++) {
        for (int j=0; j<chunk_size; j++) {
            D0[i][j] = D[i][j];
        }
    }
}

//function to retrieve the collumn and row needed from memory
void get_k_col_and_row(int n, int k, int** D0, int* k_col, int*k_row) {
    //mpi variables
    int pid, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //size variables
    int sqrt_np = sqrt(np);
    int chunk_size = n/sqrt_np;
    int offset = k/chunk_size;

    //retrieve collumns and rows needed
    for (int i=0; i<sqrt_np; i++) {
        //collumn retrieval
        if (pid==i*sqrt_np + offset) {
            for (int j=0; j<chunk_size; j++) {
                k_col[j] = D0[j][k%chunk_size];
            }
        }

        //row retrieval
        if (pid==offset*sqrt_np + i) {
            for (int j=0; j<chunk_size; j++) {
                k_row[j] = D0[k%chunk_size][j];
            }
        }
    }
}

//one to many along the rows
int rowwise_one_to_many(int* k_col, int chunk_size, int k, int row_sender_pid) {
    //mpi variables
    int pid, np;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //size variables
    int sqrt_np = sqrt(np);

    //Begin function
    int virtual_pid, iteration, mtag, reciever_pid, sender_pid;
    int distance_from_column_0 = pid % sqrt_np;
    int column_0_pid = pid - distance_from_column_0;
    int pid_row_number = (int)(pid/sqrt_np);

    int distance_sender_column_0 = row_sender_pid % sqrt_np;
    int sender_column_0 = pid - distance_from_column_0;
    int sender_row_number = (int)(pid/sqrt_np);

    //If not in sender's row, return
    if(pid_row_number == sender_row_number) {
        //Set a virtual pid to enable "negative indexing"
        if (pid < row_sender_pid){
            virtual_pid = pid + sqrt_np;
        } else {
            virtual_pid = pid;
        }

        //recieve data
        if (pid == row_sender_pid){
            iteration = 0;
        } else {
            for (int i=1; i <= log2(sqrt_np); i++){
                if (pow(2, i) + distance_sender_column_0 > virtual_pid - (pid_row_number * sqrt_np)){
                    mtag = k*sqrt_np+pid;
                    sender_pid = virtual_pid - pow(2,i-1);
                    iteration = i;
                    MPI_Recv(k_col, chunk_size, MPI_INT, sender_pid, mtag, MPI_COMM_WORLD, &status);
                    //printf("(row)%d recieves from %d \n", pid, sender_pid);
                    break;
                }
            }
        }

        //send data
        for (; iteration < log2(sqrt_np); iteration++){
            reciever_pid = pow(2, iteration) + pid;

            if (reciever_pid >= (pid_row_number + 1) * sqrt_np) {
                reciever_pid = reciever_pid - sqrt_np;
            }

            mtag = k*sqrt_np + reciever_pid;

            MPI_Send(k_col, chunk_size, MPI_INT, reciever_pid, mtag, MPI_COMM_WORLD);
            //printf("(row)%d sends to %d \n", pid, reciever_pid);
        }
    }

    return 0;
}

//one to many along collumns
int colwise_one_to_many(int* k_row, int chunk_size, int k, int column_sender_pid) {
    //mpi variables
    int pid, np;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //size variables
    int n = chunk_size*np;
    int sqrt_np = sqrt(np);

    //Begin function
    int iteration, virtual_pid, mtag, sender_pid, reciever_pid;
    int pid_column_number = pid % sqrt_np;
    int sender_column_number = column_sender_pid % sqrt_np;
    int distance_sender_row_0 = column_sender_pid - pid_column_number;

    //If not in sender's column, return
    if(pid_column_number == sender_column_number) {
        //Virtual pid to enable wrapping around
        if(pid < column_sender_pid){
            virtual_pid = pid + np;
        } else {
            virtual_pid = pid;
        }

        //recieve data
        if(pid == column_sender_pid){
            iteration = 0;
        } else {
            for (int i=1; i <= log2(sqrt_np); i++){
                if((pow(2,i) * sqrt_np - distance_sender_row_0 > virtual_pid)){
                    mtag = n*sqrt_np + k*sqrt_np + pid;
                    sender_pid = virtual_pid - pow(2,i-1) * sqrt_np;
                    MPI_Recv(k_row, chunk_size, MPI_INT, sender_pid, mtag, MPI_COMM_WORLD, &status);
                    //printf("(col)%d recieves from %d \n", pid, sender_pid);
                    iteration = i;
                    break;
                }
            }
        }

        //send data
        for (; iteration < log(sqrt_np); iteration++) {
            reciever_pid = pow(2, iteration) * sqrt_np + pid;

            if(reciever_pid > (np - sqrt_np + pid_column_number)){
                reciever_pid = reciever_pid - np;
            }

            mtag = n*sqrt_np + k*sqrt_np + reciever_pid;
            MPI_Send(k_row, chunk_size, MPI_INT, reciever_pid, mtag, MPI_COMM_WORLD);
            //printf("(col)%d sends to %d \n", pid, reciever_pid);
        }
    }
    return 0;
}

//check if a number is square
int isSquare(int n) {
    int root = sqrt(n);
    return (root*root)==n;
}

//find minimum value among 3 integers
int find_min(int a, int b, int c) {
    if (a<b && a<c)
        return a;
    if (b<c)
        return b;
    return c;
}