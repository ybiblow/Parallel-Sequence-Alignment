#pragma once

#include <stdio.h>
#include <string.h>
#include "omp.h"

#define PART  100
#define INPUT_FILE_NAME "input.txt"
#define OUTPUT_FILE_NAME "output.txt"
#define SEQ1_LENGTH 5000
#define SEQ2_LENGTH 3000
#define NUMBER_OF_WEIGHTS 4
#define SIZE_OF_COMP_MATRIX 26 * 26

char** readFromFile(char* file_name, float *weights, char* seq1, int* num_of_seq2);
void writeToFile(FILE* file, char* final_result);
char* get_Mutant_CUDA(char* sequence,int len, int m, int n);
float calc_score(char* arr, int size, float* weights);
char find_similarity(char seq1_char, char seq2_char);
char* calc_best_score_CUDA(char* seq1, char* seq2, float* comp_matrix);
void calcPortion(int* portion, int num);
void createCompMatrix(float* arr, int size, float* weights);
float findSimilarityWeight(char a, char b, float* weights);
int isIdentical(char a, char b);
int is_conservative(char a, char b);
int is_semi_conservative(char a, char b);
void CPUGetNK(int mutant_num, int seq2_len, int* n, int* k);
char* calcBestScoreOmp(float* mutantsBestScores, int* mutantsBestOffsets, int num_mutants, int seq2_len);

char* calc_similarity(char* seq1, char* seq2, int offset);
void print_char_array(char* arr, int size);
void calc_best_score(char* seq1, char* mutant, float* weights,int* best_offset, float* best_score);
void calc_num_of_occurrences(char* arr, int size, int* counter, char chr);
char* get_Mutant(char* sequence,int len, int m, int n);
