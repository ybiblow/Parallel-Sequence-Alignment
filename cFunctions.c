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
	
	*n = i;	
	*k = i + mutant_num;
}

void calcBestScoreOmp(float* mutantsBestScores, int* mutantsBestOffsets, int num_mutants, int seq2_len){
	float sscore[4] = {-100000, -100000, -100000, -100000};
	int ooffset[4];
	int mmutant[4];
	#pragma omp parallel for
	for(int i = 0; i < num_mutants; i++){
		int tid = omp_get_thread_num();
		if(mutantsBestScores[i] > sscore[tid]){
			sscore[tid] = mutantsBestScores[i];
			ooffset[tid] = mutantsBestOffsets[i];
			mmutant[tid] = i;
			
		}
	}
	//printf("score = %1.2f, offset = %d, mutant_num = %d\n", sscore[0], ooffset[0], mmutant[0]);
	//printf("score = %1.2f, offset = %d, mutant_num = %d\n", sscore[1], ooffset[1], mmutant[1]);
	//printf("score = %1.2f, offset = %d, mutant_num = %d\n", sscore[2], ooffset[2], mmutant[2]);
	//printf("score = %1.2f, offset = %d, mutant_num = %d\n", sscore[3], ooffset[3], mmutant[3]);
	float best_score = -100000;
	int best_offset;
	int mutant_num;
	for(int i = 0; i < 4; i++){
		if(sscore[i] > best_score){
			best_score = sscore[i];
			best_offset = ooffset[i];
			mutant_num = mmutant[i];
		}
	}
	int n,k;
	CPUGetNK(mutant_num + 1, seq2_len, &n, &k);
	//printf("best score = %1.2f, best offset = %d, mutant num = %d\n", best_score, best_offset, mutant_num);
	printf("mutant num: %d, MS(%d,%d), score: %1.2f, offset: %d\n", mutant_num, n, k, best_score, best_offset);
}
