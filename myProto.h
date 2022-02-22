#pragma once

#include <stdio.h>

#define PART  100
#define INPUT_FILE_NAME "input3.txt"
#define OUTPUT_FILE_NAME "output.txt"
#define SEQ1_LENGTH 5000
#define SEQ2_LENGTH 3000
#define NUMBER_OF_WEIGHTS 4

void test(int *data, int n);
int computeOnGPU(int *data, int n);
char** readFromFile(char* file_name, float *weights, char* seq1, int* num_of_seq2);
void writeToFile(FILE* file, char* file_name, int m, int n, int offset, float score);
char* get_Mutant_CUDA(char* sequence,int len, int m, int n);
void calc_best_score_CUDA(char* seq1, char* mutant, float* weights,int* best_offset, float* best_score);
char* calc_similarity_CUDA(char* seq1, char* mutant, int offset);
float calc_score(char* arr, int size, float* weights);
char find_similarity(char seq1_char, char seq2_char);
void calc_best_score_CUDA_1(char* seq1, char* mutant, float* weights,int* best_offset, float* best_score);
