#include <mpi.h> // Include MPI library for parallel processing
#include <stdio.h> // Include standard input-output library
#include <stdlib.h> // Include standard library for memory allocation
#include <limits.h> // Include library for integer limits

// Function to perform parallelized Dijkstra's algorithm
void HW2(int n, int **matrix, int *output) {
    int i, j, count, tmp, leastVal, leastPos, *done; // Declare variables
    int rank, size, local_leastVal, local_leastPos; // Declare variables for MPI communication

    // Initialize MPI rank and size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Allocate memory for array to track visited nodes
    done = (int *) calloc(n, sizeof(int));

    // Initialize output array with values from the first row of the matrix
    for (i = 0; i < n; i++) {
        done[i] = 0; // Mark all nodes as not visited
        output[i] = matrix[0][i]; // Initialize output array with values from the first row of the matrix
    }
    done[0] = 1; // Mark the first node as visited
    count = 1; // Initialize count of visited nodes

    // Loop until all nodes have been visited
    while (count < n) {
        leastVal = INT_MAX; // Initialize least value as maximum integer value
        local_leastVal = INT_MAX; // Initialize local least value as maximum integer value

        // Parallelized loop to find the local minimum value and its position
        for (i = rank; i < n; i += size) {
            tmp = output[i];
            if ((!done[i]) && (tmp < local_leastVal)) {
                local_leastVal = tmp;
                local_leastPos = i;
            }
        }

        // Reduce operation to find the global minimum value
        MPI_Allreduce(&local_leastVal, &leastVal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

        // Broadcast the position of the global minimum to all processes
        if (local_leastVal == leastVal) {
            leastPos = local_leastPos;
        }
        MPI_Bcast(&leastPos, 1, MPI_INT, rank, MPI_COMM_WORLD);

        // Mark the node with the global minimum as visited
        done[leastPos] = 1;
        count++; // Increment count of visited nodes

        // Parallelized loop to update the output array with new shortest paths
        for (i = rank; i < n; i += size) {
            if (!done[i])
                output[i] = (output[i] < leastVal + matrix[leastPos][i]) ? output[i] : leastVal + matrix[leastPos][i];
        }

        // Allgather operation to collect updated output arrays from all processes
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, output, n / size, MPI_INT, MPI_COMM_WORLD);
    }

    // Free memory allocated for the 'done' array
    free(done);
}

int main(int argc, char *argv[]) {
    int n = 4; // Size of the matrix
    int **matrix = malloc(n * sizeof(int *)); // Allocate memory for the matrix
    int *output = malloc(n * sizeof(int)); // Allocate memory for the output array
    int i, j; // Declare loop variables

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Initialize matrix and output arrays
    for (i = 0; i < n; i++) {
        matrix[i] = malloc(n * sizeof(int)); // Allocate memory for each row of the matrix
        for (j = 0; j < n; j++) {
            matrix[i][j] = (i == j) ? 0 : i + j; // Initialize matrix elements with appropriate values
        }
    }

    // Perform parallelized Dijkstra's algorithm
    HW2(n, matrix, output);

    // Finalize MPI
    MPI_Finalize();

    // Print output array
    if (output) {
        printf("Output: ");
        for (i = 0; i < n; i++) {
            printf("%d ", output[i]); // Print each element of the output array
        }
        printf("\n");
    }

    // Clean up allocated memory
    for (i = 0; i < n; i++) {
        free(matrix[i]); // Free memory allocated for each row of the matrix
    }
    free(matrix); // Free memory allocated for the matrix
    free(output); // Free memory allocated for the output array

    return 0; // Return success
}
