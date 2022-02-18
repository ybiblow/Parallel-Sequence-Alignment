#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include "myProto.h"

char* calc_similarity(char* seq1, char* seq2, int offset);
char find_similarity(char seq1_char, char seq2_char);
void print_char_array(char* arr, int size);
int is_identical(char a, char b);
int is_conservative(char a, char b);
int is_semi_conservative(char a, char b);
void calc_best_score(char* seq1, char* seq2, float* weights);
float calc_score(char* arr, int size, float* weights);
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
	/* setting variables to read from input file */
	char* file_name = (char*)FILE_NAME;
	int num_of_seq2, i, j;
	float weights[4];
	char seq1[5000];
	char** seq2 = readFromFile(file_name, &weights[0], &seq1[0], &num_of_seq2);
	
	/* calculate number of mutants and create mutants array */
	int counter = 0;
	int num_of_mutants = (strlen(seq2[0]) * (strlen(seq2[0]) - 1)) / 2;
	mutant mutants[num_of_mutants];
	for(i = 1; i < strlen(seq2[0]); i++){
		for(j = i+1; j < strlen(seq2[0]) + 1; j++){
			printf("(%d,%d) ", i, j);
			mutants[counter].m = i;
			mutants[counter].n = j;
			counter += 1;
		}
		printf("\n");
	}
		
	/* calculating the best score of seq2 */
	for(i = 0; i < num_of_seq2; i++){
		int len = strlen(seq2[i]);
		printf("i = %d\n", i);
		for(j = 0; j < num_of_mutants; j++){
			char* tmp_mutant = get_Mutant(seq2[i], len, mutants[j].m, mutants[j].n);
			printf("Mutant is = %s\n", tmp_mutant);
			//calc_best_score(seq1, tmp_mutant, weights);
		}
	}	
	
	/*
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		printf("my thread id is: %d\n", tid);
	}
	*/
	
    } else {
	//printf("my rank is %d\n", my_rank);
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

void calc_best_score(char* seq1, char* seq2, float* weights){
	int i;
	int offset = strlen(seq1) - strlen(seq2);
	int best_score = 0, tmp_score, best_offset;
	printf("offset = %d\n", offset);
	for(i = 0; i <= offset; i++){
		tmp_score = calc_score(calc_similarity(seq1, seq2, i), strlen(seq2), weights);
		if(tmp_score > best_score){
			best_score = tmp_score;
			best_offset = i;
		}
	}
	printf("Best offset is: %d\n", best_offset);

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
	int i, j = 0;
	int f_index = m - 1;
	int e_index = n - 1;
	char* mutant = (char*)malloc((len-2) * sizeof(char));
	for(i = 0; i < len; i++){
		if(i != f_index && i != e_index){
			mutant[j] = sequence[i];
			j++;
		}
	}
	printf("len = %d, j =%d\n", len, j);
	return mutant;
}
