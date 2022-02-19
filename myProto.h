#pragma once

#include <stdio.h>

#define PART  100
#define INPUT_FILE_NAME "input.txt"
#define OUTPUT_FILE_NAME "output.txt"
#define SEQ1_LENGTH 5000
#define SEQ2_LENGTH 3000
#define NUMBER_OF_WEIGHTS 4

void test(int *data, int n);
int computeOnGPU(int *data, int n);
char** readFromFile(char* file_name, float *weights, char* seq1, int* num_of_seq2);
void writeToFile(FILE* file, char* file_name, int m, int n, int offset, float score);
