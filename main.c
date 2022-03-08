#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include "myProto.h"

char* calc_similarity(char* seq1, char* seq2, int offset);
void print_char_array(char* arr, int size);
int is_identical(char a, char b);
int is_conservative(char a, char b);
int is_semi_conservative(char a, char b);
void calc_best_score(char* seq1, char* mutant, float* weights,int* best_offset, float* best_score);
void calc_num_of_occurrences(char* arr, int size, int* counter, char chr);
char* get_Mutant(char* sequence,int len, int m, int n);

typedef struct Mutant{
	int m;
	int n;
	int offset;
	float score;
}mutant;

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
	// send proc[1] portion of seq2
	int portion[2];
	calcPortion(&portion[0], num_of_seq2);
	printf("proc[0] portion = %d\n", portion[0]);
	MPI_Send(&portion[1], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	
	
	
	int conservative_matrix_size = 26 * 26;
	float conservative_matrix[conservative_matrix_size];
	createConservativeMatrix(&conservative_matrix[0], conservative_matrix_size);
	
	// calc best score, offset for a given seq2
	//void calc_best_score_CUDA(char* seq1, char* mutant, float* weights, int* best_offset, float* best_score)
	void calc_best_score_CUDA(&seq1[0], seq2[0], &weights[0]);
	
	end_time = MPI_Wtime();
	printf("Total Time: %1.4f\n", end_time - start_time);
	fclose(output_file);
	
    } else {
	int portion;
	MPI_Recv(&portion, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	printf("proc[1] portion = %d\n", portion);
	char** seq2 = (char**)malloc(portion * sizeof(char*));
	
    }
    
    MPI_Finalize();

    return 0;
}

void print_char_array(char* arr, int size){
	int i;
	for(i = 0; i < size; i++){
		printf("%c", arr[i]);
	}
	printf("\n");	
}

void calc_best_score(char* seq1, char* mutant, float* weights,int* best_offset, float* best_score){
	int i;
	int offset = strlen(seq1) - strlen(mutant);
	int tmp_score;
	*best_score = 0;
	//printf("offset = %d\n", offset);
	for(i = 0; i <= offset; i++){
		char* similarity = calc_similarity(seq1, mutant, i);
		tmp_score = calc_score(similarity, strlen(mutant), &weights[0]);
		if(tmp_score > *best_score){
			*best_score = tmp_score;
			*best_offset = i;
		}
		free(similarity);
	}
	//printf("Best offset is: %d\n", *best_offset);
	
}

float calc_score(char* arr, int size, float* weights){
	int i;
	int num_of_stars = 0, num_of_colons = 0, num_of_points = 0, num_of_spaces = 0;
	calc_num_of_occurrences(arr, size, &num_of_stars, '*');
	calc_num_of_occurrences(arr, size, &num_of_colons, ':');
	calc_num_of_occurrences(arr, size, &num_of_points, '.');
	calc_num_of_occurrences(arr, size, &num_of_spaces, ' ');
	
	return (weights[0] * num_of_stars - weights[1] * num_of_colons - weights[2] * num_of_points - weights[3] * num_of_spaces);
	
}

void calc_num_of_occurrences(char* arr, int size, int* counter, char chr){
	int i;
	for(i = 0; i < size; i++){
		if(arr[i] == chr)
			(*counter)++;
	}
}

char* calc_similarity(char* seq1, char* seq2, int offset){
	int i;
	int seq2_len = strlen(seq2);
	char* result = (char*)malloc(seq2_len * sizeof(char));
	//#pragma omp parallel for shared(seq1, seq2, offset, result) num_threads(2)
	for(i = 0; i < seq2_len; i++){
		result[i] = find_similarity(seq1[i + offset], seq2[i]);
	}
	return result;
}

char find_similarity(char seq1_char, char seq2_char){
	if(is_identical(seq1_char, seq2_char))
		return '*';
	else if(is_conservative(seq1_char, seq2_char))
		return ':';
	else if(is_semi_conservative(seq1_char, seq2_char))
		return '.';
	else
		return ' ';
}

int is_identical(char a, char b){
	if(a == b)
		return 1;
	return 0;
}

int is_conservative(char a, char b){
	int i;
	const char* conservative_groups[9] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
	for(i = 0; i < 9; i++){
		if(strchr(conservative_groups[i], a) && strchr(conservative_groups[i], b))
			return 1;
	}
	return 0;
}

int is_semi_conservative(char a, char b){
	int i;
	const char* semi_conservative_groups[11] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
	for(i = 0; i < 11; i++){
		if(strchr(semi_conservative_groups[i], a) && strchr(semi_conservative_groups[i], b))
			return 1;
	}
	return 0;
}

char* get_Mutant(char* sequence,int len, int m, int n){
	//printf("getting mutant (%d,%d) length = %d\n", m, n, len);	
	int i;
	int f_index = m - 1;
	int e_index = n - 1;
	char* mutant = (char*)malloc((len - 1) * sizeof(char));	
	#pragma omp parallel for shared(f_index, e_index, mutant)
	for(i = 0; i < len; i++){
		if(i < f_index)
			mutant[i] = sequence[i];
		else if(i > f_index && i < e_index)
			mutant[i - 1] = sequence[i];
		else if(i > e_index)
			mutant[i - 2] = sequence[i];
	}
	mutant[len - 2] = '\0';	
	return mutant;
}
