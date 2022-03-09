#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include "myProto.h"

int main(int argc, char *argv[]) {
    
	int size, my_rank;
	float weights[NUMBER_OF_WEIGHTS];
	int *data;
	int num_of_seq2;
	int counter = 0;
	int i, j, k, m, n, portion;
	int offset;
	float score;
	int best_offset;
	float best_score;
	char seq1[5000];
	int num_of_mutants;
	MPI_Status  status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
       printf("Run the example with two processes only\n");
       MPI_Abort(MPI_COMM_WORLD, __LINE__);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0) {
	printf("my rank is %d and i'm currently reading from the file\n", my_rank);
	double start_time, end_time;
	start_time = MPI_Wtime();
	/* setting variables to read from input file */
	char* input_file_name = (char*)INPUT_FILE_NAME;
	char* output_file_name = (char*)OUTPUT_FILE_NAME;
	char** seq2 = readFromFile(input_file_name, &weights[0], &seq1[0], &num_of_seq2);
	
	FILE* output_file = fopen(output_file_name, "w");
		if(output_file == NULL){
		printf("Error in opening the file\n");
		exit(1);
	}
	
	printf("num of seq2 = %d\n", num_of_seq2);	// print num of seq2
	
	// send proc[1] portion of seq2
	int portion[2];
	calcPortion(&portion[0], num_of_seq2);
	printf("proc[0] portion = %d\n", portion[0]);
	MPI_Send(&portion[1], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	
	
	
	//int comp_matrix_size = 26 * 26;
	float comp_matrix[SIZE_OF_COMP_MATRIX];
	createCompMatrix(&comp_matrix[0], SIZE_OF_COMP_MATRIX, &weights[0]);
	
	// calc best score, offset for a given seq2
	calc_best_score_CUDA(&seq1[0], seq2[0], comp_matrix);
	
	end_time = MPI_Wtime();
	printf("Total Time: %1.4f\n", end_time - start_time);
	fclose(output_file);
	
    } else {
	int portion;
	MPI_Recv(&portion, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	//printf("proc[1] portion = %d\n", portion);
	//char** seq2 = (char**)malloc(portion * sizeof(char*));
	
    }
    
    MPI_Finalize();

    return 0;
}
