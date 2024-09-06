/*
Group 5
Brandon Bratcher    R#11520450
James Gibson        R#11671129
Richard Pearson     R#11786900
Dhruv Maniar        R#11713343
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

//constants
#define c1 1.23456
#define c2 6.54321
#define n 16

//function declarations
void hw5(double*,double*);
void many_to_one(int,double*,double*);
void default_x(int,double*);
void default_f(int,double*);
void sum_arrays(int,double*,double*,double*);
int is_divis_by_twice(int,int);
int isPowerOf2(int);


int main(int argc, char **argv) {
    //variable declaration
    int pid, np;
    //MPI variable initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //check if number of processes is valid
    //processes must be a power of two
    if (!isPowerOf2(np)) {
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
        MPI_Finalize();
        return -1;
    }
    //n must be divisible by 2p
    if (!is_divis_by_twice(n,np)) {
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
        MPI_Finalize();
        return -1;
    }

    //create the position array
    double *x = (double*)malloc(n * sizeof(double));
    default_x(n, x);
    
    //create the force array
    double *f = (double*)malloc(n * sizeof(double));
    default_f(n,f);

    //calculate the force on every particle
    hw5(x,f);
    //combine all calculated forces into one process
    double *temp_f = (double*)malloc(n * sizeof(double));
    many_to_one(n,f,temp_f);

    //print output
    if (pid == 0){
        for (int i=0; i<n; i++){
            //print force
            printf("%3.3f ", f[i]);
            
            //break line
            if(i%10 == 0 && i != 0){
                printf("\n");
            }
        }
    }

    //exit
    MPI_Finalize();
    return 0;
}

//calculates the force applied on each particle from each other particle
// (parrellized so that each process has equivalent workload)
void hw5(double *x, double *f) {
    //variable declaration
    int pid, np;
    //MPI variable initialization
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int gap_size = n/(2*np);
    
    //perform program on an earlier section of the array
    for (int i=(pid*gap_size)+1; i<=((pid+1)*gap_size); i++) {
        for (int j=0; j<i; j++) {
            double diff = x[i] - x[j];
            double temp = 1.0/(diff*diff*diff);
            double tmp = c1*temp*temp - c2*temp;
            f[i] += tmp;
            f[j] -= tmp;
        }
    }
    
    //perform program on a later section of the array
    for (int i=n-((pid+1)*gap_size)+1; i<=n-(pid*gap_size); i++) {
        for (int j=0; j<i; j++) {
            double diff = x[i] - x[j];
            double temp = 1.0/(diff*diff*diff);
            double tmp = c1*temp*temp - c2*temp;
            f[i] += tmp;
            f[j] -= tmp;
        }
    }
}

//communication process to send all info back to a single process
void many_to_one(int length, double* data_send, double* data_recv) {
    //variable declaration
    int pid, np, mtag, mtag1, mtag2, reciever_pid, sender_pid;
    MPI_Status status;
    //MPI variable initialization
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);


    for (int i=0; i<np/2;i++) {
        mtag = i;
        //if the process is odd, send. Else if the process is even, recieve.
        //Creates this pattern O X O X O X O X
        if (pid == 2*i+1) {
            reciever_pid = pid-1;
            MPI_Send(data_send, length, MPI_DOUBLE, reciever_pid, mtag, MPI_COMM_WORLD);
        } else if (pid == 2*i) {
            sender_pid = pid+1;
            MPI_Recv(data_recv, length, MPI_DOUBLE, sender_pid, mtag, MPI_COMM_WORLD, &status);
            //after recieving, add values to sending array
            sum_arrays(length,data_send,data_send,data_recv);
        }
    }
    mtag1 = np/2;
    mtag2 = np/2;

    //For the number of send/recieves required
    for (int j=1; j<log2(np); j++) {
        //For the number of send calls in this iteration
       for (int i=1;i<np/pow(2,j);i+=2) {
            /**
             * X indicate process that have sent data
             * O X O X O X O X
             * O X X X O X X X
             * O X X X X X X X < This line indicates that all processes have sent to process 0
             **/
            if (pid == pow(2,j)*i) {
                reciever_pid = pow(2,j)*(i-1);
                MPI_Send(data_send, length, MPI_DOUBLE, reciever_pid, mtag1, MPI_COMM_WORLD);
            } 

            mtag1++;
        }

        //This is the recieve loop, following the same process as above
        for (int i=0;i<np/pow(2,j);i+=2) {
            if(pid == pow(2,j)*i) {
                sender_pid = pow(2,j)*(i+1);
                MPI_Recv(data_recv, length, MPI_DOUBLE, sender_pid, mtag2, MPI_COMM_WORLD, &status);
                //after recieving, add values to sending array
                sum_arrays(length,data_send,data_send,data_recv);
            }

            mtag2++;
        }
    }
}

//generate default positioon array
void default_x(int size, double *x) {
    for (int i=0; i<size; i++) {
        x[i] = i;
    }
}

//generate default force array
void default_f(int size, double *f) {
    for (int i=0; i<size; i++) {
        f[i] = 0;
    }
}

//add two arrays together
void sum_arrays(int length, double *arr1, double *arr2, double *arr3) {
    for (int i=0; i<length; i++) {
        arr1[i] = arr2[i] + arr3[i];
    }
}

//check if a is divisible by 2b
int is_divis_by_twice(int a, int b) {
    float division = a/(2*b);
    float floored_division = floor(division);
    return division == floored_division;
}

//check if a number is a power of 2
int isPowerOf2(int np) {
    float log2np = log2(np);
    float floor_log2np = floor(log2np);
    return log2np==floor_log2np;
}
