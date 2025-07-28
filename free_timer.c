#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ALLOC_SIZE 5 // Each allocation is 5 bytes

int main(int argc, char *argv[]) {
    int rank, size;
    
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Check if the correct number of arguments are provided
    if (argc != 2) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <num_allocs>\n", argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // Parse the number of allocations from the command-line argument
    long NUM_ALLOCS = atol(argv[1]);
    
    // Allocate an array of pointers
    void **pointers = (void **)malloc(NUM_ALLOCS * sizeof(void *));
    if (pointers == NULL) {
        fprintf(stderr, "Process %d: Failed to allocate memory for pointer array\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Perform malloc calls
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    for (long i = 0; i < NUM_ALLOCS; i++) {
        pointers[i] = malloc(ALLOC_SIZE);
        if (pointers[i] == NULL) {
            fprintf(stderr, "Process %d: Failed to malloc at iteration %ld\n", rank, i);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    double end_time = MPI_Wtime();
    double malloc_time = end_time - start_time;

    // Synchronize processes before starting the timer for freeing memory
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // Free all allocated memory
    for (long i = 0; i < NUM_ALLOCS; i++) {
        free(pointers[i]);
    }

    // End timer
    end_time = MPI_Wtime();
    double free_time = end_time - start_time;

    // Free the pointer array
    free(pointers);

    // Print time taken to free memory
    printf("Process %d: Time taken to malloc %ld allocations = %f seconds\n", rank, NUM_ALLOCS, malloc_time);
    printf("Process %d: Time taken to free %ld allocations = %f seconds\n", rank, NUM_ALLOCS, free_time);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}