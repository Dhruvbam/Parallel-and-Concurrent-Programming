#include <stdio.h>
#include <math.h>
#include <mpi.h>

// Custom barrier function to synchronize MPI processes
void mybarrier(MPI_Comm comm) {
    int rank, size;
    // Get the rank (ID) of the current process
    MPI_Comm_rank(comm, &rank);
    // Get the total number of processes
    MPI_Comm_size(comm, &size);

    // Calculate the number of steps required for the barrier
    // This is based on the number of processes, and it's the ceiling of log2(size)
    int steps = (int)ceil(log2(size));
    MPI_Status status;

    // First phase: Each process sends a message to its partner in each step
    for (int i = 0; i < steps; i++) {
        // Calculate the partner's rank for this step
        // The partner is determined by XORing the current rank with 2^i
        int partner = rank ^ (1 << i);
        // Check if the partner is within the valid range of processes
        if (partner < size) {
            // Send a message to the partner
            MPI_Send(&rank, 1, MPI_INT, partner, 0, comm);
            // Receive a message from the partner (synchronization point)
            MPI_Recv(NULL, 0, MPI_BYTE, partner, 0, comm, &status);
        }
    }

    // Second phase: Each process sends a confirmation message back to its partners in reverse order
    for (int i = steps - 1; i >= 0; i--) {
        // Calculate the partner's rank for this step (same as in the first phase)
        int partner = rank ^ (1 << i);
        // Check if the partner is within the valid range of processes
        if (partner < size) {
            // Send a confirmation message to the partner
            MPI_Send(NULL, 0, MPI_BYTE, partner, 0, comm);
            // Receive a confirmation message from the partner (synchronization point)
            MPI_Recv(&partner, 1, MPI_INT, partner, 0, comm, &status);
        }
    }
}

int main(int argc, char *argv[]) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Call the custom barrier function to synchronize all processes
    mybarrier(MPI_COMM_WORLD);

    // Finalize the MPI environment
    MPI_Finalize();
    return 0;
}
