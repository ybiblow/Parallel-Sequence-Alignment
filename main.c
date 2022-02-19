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
void calc_best_score(char* seq1, char* seq2, float* weights,int* best_offset, float* best_score);
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
	int num_of_seq2, i, j, k;
	float weights[4];
	char seq1[5000];
	char** seq2 = readFromFile(file_name, &weights[0], &seq1[0], &num_of_seq2);
	
	/* calculate number of mutants and create mutants array */
	int counter = 0;
	int num_of_mutants = (strlen(seq2[0]) * (strlen(seq2[0]) - 1)) / 2;
	printf("number of mutants = %d\n", num_of_mutants);
	//mutant* mutants[num_of_seq2];
	mutant** mutants = (mutant**)malloc(num_of_seq2 * sizeof(mutant*));
	
	for(i = 0; i < num_of_seq2; i++){
		num_of_mutants = (strlen(seq2[i]) * (strlen(seq2[i]) - 1)) / 2;
		printf("SEQ2 num = %d, num of mutants = %d\n", i, num_of_mutants);
		mutants[i] = (mutant*)malloc(num_of_mutants * sizeof(mutant));
		counter = 0;
		for(j = 1; j < strlen(seq2[i]); j++){
			for(k = j + 1; k < strlen(seq2[i]) + 1; k++){
				printf("(%d,%d) ", j, k);
				int offset;
				float score;
				mutants[i][counter].m = j;
				mutants[i][counter].n = k;
	
				char* tmp_mutant = get_Mutant(seq2[i], strlen(seq2[i]), mutants[i][counter].m, mutants[i][counter].n);
				printf("SEQ2_%d, mutant_%d = %s\n", i, counter, tmp_mutant);
	
				/* calculating the best score of each mutant  in seq2 */
				calc_best_score(seq1, tmp_mutant, weights, &offset, &score);
				mutants[i][counter].score = score;
				mutants[i][counter].offset = offset;
				counter++;
				
			}
			printf("\n");
		}
	}
	
	for(i = 0; i < num_of_seq2; i++){
		printf("SEQ2_%d results:\n", i);
		int len = strlen(seq2[i]);
		num_of_mutants = len * (len -1) / 2;
		for(j = 0; j < num_of_mutants; j++){
			printf("mutant (%d,%d) offset = %d score %f\n", mutants[i][j].m, mutants[i][j].n, mutants[i][j].offset, mutants[i][j].score);
		}
	}
	
	for(i = 0; i < num_of_seq2; i++){
		free(mutants[i]);
	}
	
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

void calc_best_score(char* seq1, char* seq2, float* weights,int* best_offset, float* best_score){
	int i;
	int offset = strlen(seq1) - strlen(seq2);
	int tmp_score;
	*best_score = 0;
	//printf("offset = %d\n", offset);
	for(i = 0; i <= offset; i++){
		tmp_score = calc_score(calc_similarity(seq1, seq2, i), strlen(seq2), &weights[0]);
		if(tmp_score > *best_score){
			*best_score = tmp_score;
			*best_offset = i;
		}
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
	
	int i, j = 0;
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
	
	//printf("mutant length = %d\n", (int)strlen(mutant));
	//printf("mutant after parallel = %s\n", mutant);
	
	/*
	for(i = 0; i < len; i++){
		if(i != f_index && i != e_index){
			mutant[j] = sequence[i];
			j++;
		}
	}
	mutant[j] = '\0';
	printf("mutant after sequential = %s\n", mutant);*/
	
	//printf("len = %d, j =%d, mutant len = %d\n", len, j, (int)strlen(mutant));
	
	return mutant;
}
