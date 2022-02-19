#pragma once

#define PART  100
#define FILE_NAME "input.txt"
#define SEQ1_LENGTH 5000
#define SEQ2_LENGTH 3000
#define NUMBER_OF_WEIGHTS 4

void test(int *data, int n);
int computeOnGPU(int *data, int n);
char** readFromFile(char* file_name, float *weights, char* seq1, int* num_of_seq2);
