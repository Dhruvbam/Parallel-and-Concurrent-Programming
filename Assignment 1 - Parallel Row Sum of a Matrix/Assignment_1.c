#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> // Include the MPI library for parallel programming

#define generate_data(i, j) ((i) + (j) * (j)) // Macro to generate data for the matrix

int main(int argc, char **argv)
{
    int i, j, pid, np, mtag;
    int data[100][100], row_sum[50]; // Declare arrays for data and row sums
    MPI_Status status; // Status variable for MPI communication
    MPI_Request req_s, req_r; // Request variables for non-blocking communication
    int flag = 0; // Flag for checking completion of MPI_Irecv

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // Get the rank (ID) of the current process
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    // Get the total number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // Process 0 generates and sends data, and computes row sums for its half
    if (pid == 0)
    {
        // Generate data for rows 0-49 and send to pid 1
        for (i = 0; i < 50; i++)
            for (j = 0; j < 100; j++)
                data[i][j] = generate_data(i, j); // Fill the data array

        mtag = 1; // Message tag for communication
        // Start a non-blocking send of the data to process 1
        MPI_Isend(data, 5000, MPI_INT, 1, mtag, MPI_COMM_WORLD, &req_s);

        // Generate data and compute row sums for rows 50-99
        for (i = 50; i < 100; i++)
        {
            row_sum[i - 50] = 0;
            for (j = 0; j < 100; j++)
            {
                data[i][j] = generate_data(i, j);
                row_sum[i - 50] += data[i][j]; // Compute row sum
            }
        }

        // Wait for the send operation to complete
        MPI_Wait(&req_s, &status);

        // Receive computed row sums from pid 1
        mtag = 2;
        MPI_Recv(row_sum + 50, 50, MPI_INT, 1, mtag, MPI_COMM_WORLD, &status);

        // Print all row sums
        for (i = 0; i < 100; i++)
        {
            printf("%d ", row_sum[i]);
            if ((i + 1) % 10 == 0)
                printf("\n");
        }
    }
    else // pid == 1 receives data and computes row sums for its half
    {
        mtag = 1;
        // Start a non-blocking receive of the data from process 0
        MPI_Irecv(data, 5000, MPI_INT, 0, mtag, MPI_COMM_WORLD, &req_r);

        // Compute row sums for the first 25 rows while waiting for data
        for (i = 0; i < 25; i++)
        {
            row_sum[i] = 0;
            for (j = 0; j < 100; j++)
                row_sum[i] += data[i][j]; // Compute row sum
        }

        // Check if the receive operation has completed
        MPI_Test(&req_r, &flag, &status);

        // Compute row sums for the remaining rows as data becomes available
        for (i = 25; i < 50; i++)
        {
            // Wait for data to be available if not already received
            if (!flag)
                MPI_Wait(&req_r, &status);

            row_sum[i] = 0;
            for (j = 0; j < 100; j++)
                row_sum[i] += data[i][j]; // Compute row sum
        }

        // Send computed row sums to pid 0
        mtag = 2;
        MPI_Send(row_sum, 50, MPI_INT, 0, mtag, MPI_COMM_WORLD);
    }

    // Finalize the MPI environment
    MPI_Finalize();
    return 0;
}
