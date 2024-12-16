#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_ALLOCS 55213440
#define ALLOC_SIZE 5 // Each allocation is 5 bytes (5 * NUM_ALLOCS = 280 MB)
#define GAP_SIZE 63 * 5 
int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Allocate an array of pointers
    void **pointers = (void **)malloc(NUM_ALLOCS * sizeof(void *));
    if (pointers == NULL) {
        fprintf(stderr, "Process %d: Failed to allocate memory for pointer array\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Allocate an array of gaps
    void **gaps = (void **)malloc(NUM_ALLOCS * sizeof(void *));
    if (gaps == NULL) {
        fprintf(stderr, "Process %d: Failed to allocate memory for gaps array\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Perform malloc calls with gaps
    for (long i = 0; i < NUM_ALLOCS; i++) {
        pointers[i] = malloc(ALLOC_SIZE);
        if (pointers[i] == NULL) {
            fprintf(stderr, "Process %d: Failed to malloc at iteration %ld\n", rank, i);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        gaps[i] = malloc(GAP_SIZE);
        if (gaps[i] == NULL) {
            fprintf(stderr, "Process %d: Failed to malloc at iteration %ld\n", rank, i);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    //free gaps first
    for (long i = 0; i < NUM_ALLOCS; i++) {
        free(gaps[i]);
    }
    free(gaps);

    // Synchronize processes before starting the timer
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // Free all allocated memory
    for (long i = 0; i < NUM_ALLOCS; i++) {
        free(pointers[i]);
    }

    // End timer
    double end_time = MPI_Wtime();
    double free_time = end_time - start_time;

    // Free the pointer array
    free(pointers);

    // Print time taken to free memory
    printf("Process %d: Time taken to free %ld allocations = %f seconds\n", rank, NUM_ALLOCS, free_time);

    MPI_Finalize();
    return 0;
}