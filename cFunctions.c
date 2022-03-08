#include "myProto.h"
#include <stdio.h>
#include <stdlib.h>

void test(int *data, int n) {
    int i;
    for (i = 0;   i < n;   i++) {
        if (data[i] != i + 1) {
           printf("Wrong Calculations - Failure of the test at data[%d]\n", i);
           return;
    	}
    }
    printf("The test passed successfully\n"); 
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

void writeToFile(FILE* file, char* file_name, int m, int n, int offset, float score){
	if(file == NULL){
		printf("writeToFile() - file is NULL\n");
		exit(1);
	}
		fprintf(file, "MS(%d,%d), Offset = %d, Score = %1.2f\n",m ,n, offset, score);
}

void calcPortion(int* portion, int num){
	portion[0] = num /2;
	portion[1] = (num % 2 != 0) ? portion[1] = num / 2 + 1 : portion[1] = num / 2;
}

void createConservativeMatrix(float* arr, int size){
	for(int i = 0; i < 26; i++){
		for(int j = 0; j < 26; j++){
			conservative_matrix[i * 26 + j] = findSimilarityWeight(i + 'A', j + 'A', &weights[0]);
			//printf("%1.2f ", conservative_matrix[i * 26 + j]);
		}
		//printf("\n");
	}
}
