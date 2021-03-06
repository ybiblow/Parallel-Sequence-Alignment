#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include "myProto.h"

int main(int argc, char *argv[]) {
    
	int size, my_rank;
	float weights[NUMBER_OF_WEIGHTS];
	int num_of_seq2;
	char seq1[SEQ1_LENGTH];
	int num_of_mutants;
	MPI_Status  status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
       printf("Run the example with two processes only\n");
       MPI_Abort(MPI_COMM_WORLD, __LINE__);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    /*
    	proc[0] is the master and does most of the fundumental tasks like reading, printing from/to a file and sending,receiving the relevant information from/to proc[1].

    	proc[1] receives information from proc[0] and works on a given number of seq2's, when finished proc[1] send the results to proc[0].
    	
    	both procs use the same function to calc the best score, offset and (n,k) for a given seq2, that function uses CUDA to calc the results for each mutant and we use
    	OMP to find the best result out of the many results for each mutant. it's important to note that the results for the scores and offsets for those results are store in 2 vectors with the
    	length of the number of mutants for a given sequence2.
    	
    	both procs use the compare matrix for getting the score of two chars (one from seq1 and another from seq2).
    	
    */
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
	
	// send proc1[1] seq1
	int seq1_len = strlen(seq1) + 1;
	MPI_Send(&seq1_len, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	MPI_Send(&seq1[0], seq1_len, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
	
	// send proc[1] weights
	MPI_Send(&weights, 4, MPI_INT, 1, 0, MPI_COMM_WORLD);
	
	// calculating proc[0] and proc[1] portions of the number of seq2's and sending proc[1] its portion so it can allocate mem for receiving seq2's
	int portion[2];
	calcPortion(&portion[0], num_of_seq2);
	printf("proc[0] portion = %d\n", portion[0]);	
	MPI_Send(&portion[1], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	
	// creating comp_matrix - a matrix for comparing 2 chars and finding out if they are equal of in a conservative group etc
	float comp_matrix[SIZE_OF_COMP_MATRIX];
	createCompMatrix(&comp_matrix[0], SIZE_OF_COMP_MATRIX, &weights[0]);
	
	// sending proc[1] his portions of seq's
	for(int i = portion[0]; i < num_of_seq2; i++){
		int len = strlen(seq2[i]) + 1;
		MPI_Send(&len, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Send(seq2[i], len, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
	}
	
	// for each seq2 from proc[0] portion calculating best score, offset and (N,K)
	for(int i = 0; i < portion[0]; i++){
		//returns a string that is to be printed in the output file
		char* final_result = calc_best_score_CUDA(&seq1[0], seq2[i], comp_matrix);
		writeToFile(output_file, final_result);
		free(final_result);
	}
	// receiving from proc[1] it's portion of results for each of the seq2's
	for(int i = 0; i < portion[1]; i++){
		int tmp_length;
		char* tmp;
		MPI_Recv(&tmp_length, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
		tmp = (char*) malloc(tmp_length * sizeof(char));
		MPI_Recv(tmp, tmp_length, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);
		writeToFile(output_file, tmp);
		free(tmp);
	}
	
	end_time = MPI_Wtime();
	printf("Total Time: %1.4f\n", end_time - start_time);
	fclose(output_file);
	
    } else {
	int portion, seq1_length;
	char* seq1;
	// receive seq1 from proc[0]
	MPI_Recv(&seq1_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	seq1 = (char*)malloc(seq1_length * sizeof(char));
	MPI_Recv(seq1, seq1_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);

	// recv weights array from proc[0]
	MPI_Recv(&weights, 4, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);

	// recv the size of the portion of seq2's from proc[0]
	MPI_Recv(&portion, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	printf("proc[1] portion = %d\n", portion);
	char** seq2 = (char**)malloc(portion * sizeof(char*));
	
	// recv portion of seq2's from proc[0]
	for(int i = 0; i < portion; i++){
		int tmp_length;
		MPI_Recv(&tmp_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		char* tmp = (char*)malloc(tmp_length * sizeof(char));
		MPI_Recv(tmp, tmp_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
		seq2[i] = tmp;
	}
	
	float comp_matrix[SIZE_OF_COMP_MATRIX];
	createCompMatrix(&comp_matrix[0], SIZE_OF_COMP_MATRIX, &weights[0]);
	// for each seq2 calculating the best score, offset and (N,K)
	for(int i = 0; i < portion; i++){
		// final_result will contain a string that describes the best score, offset etc.
		char* final_result = calc_best_score_CUDA(&seq1[0], seq2[i], comp_matrix);
		int final_result_length = strlen(final_result) + 1;
		MPI_Send(&final_result_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&final_result[0], final_result_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
		free(final_result);
	}
	
	for(int i = 0; i < portion; i++){
		free(seq2[i]);
	}
	free(seq2);
    }
    
    MPI_Finalize();

    return 0;
}
