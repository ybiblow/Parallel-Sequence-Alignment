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
	printf("File name is: %s\n", file_name);
	FILE* file = fopen(file_name, "r");
		if(file == NULL){
			printf("Error in opening the file\n");
			exit(1);
		}
	for(i = 0; i < 4; i++)
			fscanf(file,"%f", &weights[i]);
	fscanf(file, "%s", &seq1[0]);
	printf("SEQ1 is: %s\n", &seq1[0]);
	fscanf(file, "%d", num_of_seq2);
	printf("Number of Sequences is: %d\n", *num_of_seq2);
	seq2 = (char**)malloc(*num_of_seq2 * sizeof(char*));	
	for(i = 0; i < *num_of_seq2; i++)
		seq2[i] = (char*)malloc(SEQ2_LENGTH* sizeof(char));
	for(i = 0; i < *num_of_seq2; i++)
		fscanf(file, "%s", seq2[i]);
	for(i = 0; i < *num_of_seq2; i++)
		printf("SEQ2 num %d is: %s\n", i, seq2[i]);
	
	fclose(file);	
	return seq2;
}
