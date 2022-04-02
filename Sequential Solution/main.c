#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INPUT_FILE_NAME "input3.txt"
#define OUTPUT_FILE_NAME "output.txt"
#define SEQ1_LENGTH 5000
#define SEQ2_LENGTH 3000
#define NUMBER_OF_WEIGHTS 4
#define SIZE_OF_COMP_MATRIX 26 * 26

char** readFromFile(char* file_name, float *weights, char* seq1, int* num_of_seq2);
void createCompMatrix(float* arr, int size, float* weights);
float findSimilarityWeight(char a, char b, float* weights);
int isIdentical(char a, char b);
int is_conservative(char a, char b);
int is_semi_conservative(char a, char b);
void CPUGetNK(int mutant_num, int seq2_len, int* n, int* k);

int main(int argc, char *argv[]){
	
	float weights[NUMBER_OF_WEIGHTS];
	int num_of_seq2;
	char seq1[SEQ1_LENGTH];
	int num_of_mutants;
	char* input_file_name = (char*)INPUT_FILE_NAME;
	char* output_file_name = (char*)OUTPUT_FILE_NAME;
	FILE* output_file = fopen(output_file_name, "w");
		if(output_file == NULL){
		printf("Error in opening the file\n");
		exit(1);
	}
	char** seq2 = readFromFile(input_file_name, &weights[0], &seq1[0], &num_of_seq2);
	float comp_matrix[SIZE_OF_COMP_MATRIX];
	createCompMatrix(&comp_matrix[0], SIZE_OF_COMP_MATRIX, &weights[0]);
	int i, j, s, comp_idx, n, k;
	int seq1_len = strlen(seq1);
	int maxOffset;
	int counter;
	float score;
	float best_score;
	int best_offset;
	int best_n, best_k;
	for(i = 0; i < num_of_seq2; i++){
		int seq2_len = strlen(seq2[i]);
		char* current_seq2 = seq2[i];
		maxOffset = seq1_len - (seq2_len - 2) + 1;
		num_of_mutants = seq2_len * (seq2_len - 1) / 2;
		best_score = -100000;
		for(j = 0; j < num_of_mutants; j++){
			CPUGetNK(j + 1, seq2_len, &n, &k);
			for(s = 0; s < maxOffset; s++){
				score = 0;
				counter = 0;
				for(comp_idx = 0; comp_idx < seq2_len - 2; comp_idx++, counter++){
					if(counter == n - 1)
						counter++;
					if(counter == k - 1)
						counter++;
					score += comp_matrix[(seq1[s + comp_idx] - 'A') * 26 + (current_seq2[counter] - 'A')];
					float tmp = comp_matrix[(seq1[s + comp_idx] - 'A') * 26 + (current_seq2[counter] - 'A')];
					/*
					if(n == 1 && k ==2 && s == 2){
						//printf("im here\n");
						printf("comparint between %c and %c score is: %1.4f\n", seq1[s + comp_idx], current_seq2[counter], tmp);
					}
					*/
				}
				
				//if(n == 1 && k ==2)
					//printf("im here\n");
				
				//printf("for (%d,%d) i have found score: %1.4f and offset is: %d\n", n, k, score, s);
				if(score > best_score){
					best_score = score;
					best_offset = s;
					best_n = n;
					best_k = k;
				}				
			}
		}
		printf("MS(%d,%d)	best score: %1.4f	best offset: %d\n", best_n, best_k, best_score, best_offset);	
	}
	
	
	
	printf("Hello World!\n");
	return 0;
}

char** readFromFile(char* file_name, float *weights, char* seq1, int* num_of_seq2){
	int i;
	char** seq2;
	FILE* file = fopen(file_name, "r");
		if(file == NULL){
			printf("Error in opening the file\n");
			exit(1);
		}
	for(i = 0; i < 4; i++)
			fscanf(file,"%f", &weights[i]);
	fscanf(file, "%s", &seq1[0]);
	fscanf(file, "%d", num_of_seq2);
	seq2 = (char**)malloc(*num_of_seq2 * sizeof(char*));	
	for(i = 0; i < *num_of_seq2; i++)
		seq2[i] = (char*)malloc(SEQ2_LENGTH* sizeof(char));
	for(i = 0; i < *num_of_seq2; i++)
		fscanf(file, "%s", seq2[i]);
	
	fclose(file);	
	return seq2;
}

void createCompMatrix(float* arr, int size, float* weights){
	//printf("Comp Matrix is:\n");
	for(int i = 0; i < 26; i++){
		for(int j = 0; j < 26; j++){
			arr[i * 26 + j] = findSimilarityWeight(i + 'A', j + 'A', &weights[0]);
			//printf("%1.2f ", arr[i * 26 + j]);
		}
		//printf("\n");
	}
}

float findSimilarityWeight(char a, char b, float* weights){
    if (isIdentical(a, b))
        return weights[0];
    else if(is_conservative(a, b))
        return -weights[1];
    else if(is_semi_conservative(a, b))
        return -weights[2];
    else
        return -weights[3];
}

int isIdentical(char a, char b){
    if(a == b)
        return 1;
    return 0;
}

int is_conservative(char a, char b){
    const char* conservative_groups[9] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
	for(int i = 0; i < 9; i++){
		if(strchr((char*)conservative_groups[i], a) && strchr((char*)conservative_groups[i], b))
			return 1;
	}
	return 0;
}

int is_semi_conservative(char a, char b){
    const char* semi_conservative_groups[11] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
	for(int i = 0; i < 11; i++){
		if(strchr((char*)semi_conservative_groups[i], a) && strchr((char*)semi_conservative_groups[i], b))
			return 1;
	}
	return 0;
}

void CPUGetNK(int mutant_num, int seq2_len, int* n, int* k){
	int i;
	int num_of_mutants_in_row = seq2_len;

	for(i = 1; i < seq2_len; i++){
		if(mutant_num - (num_of_mutants_in_row - 1) > 0){
		    mutant_num -= (num_of_mutants_in_row - 1);
		    num_of_mutants_in_row--;
		}else{
		    break;
		}
	}
	//printf("(%d,%d)\n", i, i + mutant_num);
	*n = i;	
	*k = i + mutant_num;
}
