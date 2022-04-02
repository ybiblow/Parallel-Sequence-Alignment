#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "myProto.h"

#define CUDA_ERR_CHECK(err,msg) {if (err != cudaSuccess) {fprintf(stderr, msg " - %s\n", cudaGetErrorString(err));exit(EXIT_FAILURE);}}

#define CUDA_MEM_INIT(d_pointer, size, type) {cudaError_t err = cudaSuccess;\
					size_t arrSize = size * sizeof(type);\
				err = cudaMalloc((void**)&d_pointer, arrSize);\
			CUDA_ERR_CHECK(err, "Failed to allocate device memory");}
	
#define CUDA_MEM_INIT_COPY(dest, src, size, type) {\
			cudaError_t err = cudaSuccess;\
		size_t  arrSize = size * sizeof(type);\
		err = cudaMalloc((void**)&dest, arrSize);\
CUDA_ERR_CHECK(err, "Failed to allocate device memory");\
err = cudaMemcpy(dest, src, arrSize, cudaMemcpyHostToDevice);\
CUDA_ERR_CHECK(err, "Failed to copy data from host to device"); }

__device__ int strchr_CUDA(char* arr, char a){
	for(;;arr++){
		if(*arr == a)
			return 1;
		if(*arr == '\0')
			return 0;
	}
	
}

__device__ int is_identical_CUDA(char a, char b){
	if(a == b)
		return 1;
	return 0;
}

__device__ int is_conservative_CUDA(char a, char b){
	int i;
	const char* conservative_groups[9] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
	for(i = 0; i < 9; i++){
		if(strchr_CUDA((char*)conservative_groups[i], a) && strchr_CUDA((char*)conservative_groups[i], b))
			return 1;
	}
	return 0;
}

__device__ int is_semi_conservative_CUDA(char a, char b){
	int i;
	const char* semi_conservative_groups[11] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
	for(i = 0; i < 11; i++){
		if(strchr_CUDA((char*)semi_conservative_groups[i], a) && strchr_CUDA((char*)semi_conservative_groups[i], b))
			return 1;
	}
	return 0;
}

__device__ char find_similarity_Kernel(char a, char b) {
	if(is_identical_CUDA(a, b))
		return '*';
	else if(is_conservative_CUDA(a, b))
		return ':';
	else if(is_semi_conservative_CUDA(a, b))
		return '.';
	else
		return ' ';
}

__global__ void calc_similarity_CUDA_Kernel(char *d_seq1, char* d_mutant, int d_mutant_len, char* d_result, int offset) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	// Increment the proper value of the arrray according to thread ID 
	if (i < d_mutant_len){
		d_result[i] = find_similarity_Kernel(d_seq1[i + offset], d_mutant[i]);
	}
}

__device__	void calc_score_cuda(char* seq1, char* mutant, float* score, float* weights, int mutant_len){
	int i, stars = 0, colons = 0, dots = 0, spaces = 0;
	char* similarity_arr = (char*)malloc(mutant_len * sizeof(char));
	for(i = 0; i < mutant_len; i++){
		similarity_arr[i] = find_similarity_Kernel(seq1[i], mutant[i]);
		if(similarity_arr[i] == '*')
			stars++;
		else if(similarity_arr[i] == ':')
			colons++;
		else if(similarity_arr[i] == '.')
			dots++;
		else if(similarity_arr[i] == ' ')
			spaces++;
	}
	free(similarity_arr);
	*score = weights[0] * stars - weights[1] * colons - weights[2] * dots - weights[3] * spaces;
}

__global__	void calcBestScore(char* d_seq1, char* d_mutant, float* d_score_arr, float* weights, int offset, int mutant_len) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	// Increment the proper value of the arrray according to thread ID 
	if (i <= offset){
		calc_score_cuda(&d_seq1[i], d_mutant, &d_score_arr[i], weights, mutant_len);
	}
}

__global__ void get_Mutant_CUDA_Kernel(char *sequence, char* d_tmp_mutant, int len, int m, int n) {
		int i = blockDim.x * blockIdx.x + threadIdx.x;
		//char* tmp_mutant = (char*)malloc((len-2)*sizeof(char));
		if (i < len){
			if(i < (m - 1))
				d_tmp_mutant[i] = sequence[i];
			else if(i > (m - 1) && i < (n - 1))
				d_tmp_mutant[i - 1] = sequence[i];
			else if(i > (n - 1))
				d_tmp_mutant[i - 2] = sequence[i];
		}
		if(i == len)
			d_tmp_mutant[i - 2] = '\0';
}

char* get_Mutant_CUDA(char* sequence,int len, int m, int n){

	char* mutant = (char*)malloc((len - 2) * sizeof(char));
	// Error code to check return values for CUDA calls
	cudaError_t err = cudaSuccess;

	size_t size = len * sizeof(char);
  

	// Allocate memory on GPU to copy the data from the host
	char* d_tmp_mutant;
	err = cudaMalloc((void **)&d_tmp_mutant, (len - 2) * sizeof(char));      
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
    
	// Allocate memory on GPU to copy the data from the host
	char* d_sequence;
	err = cudaMalloc((void **)&d_sequence, size);      
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	// Copy data from host to the GPU memory
	err = cudaMemcpy(d_sequence, sequence, size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}


	// Launch the Kernel
	int threadsPerBlock = 256;
	int blocksPerGrid =(len + threadsPerBlock - 1) / threadsPerBlock;
	get_Mutant_CUDA_Kernel<<<blocksPerGrid, threadsPerBlock>>>(d_sequence, d_tmp_mutant, len, m, n);
	err = cudaGetLastError();
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to launch vectorAdd kernel -  %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	// Copy the  result from GPU to the host memory.
	err = cudaMemcpy(mutant, d_tmp_mutant, (len - 2) * sizeof(char), cudaMemcpyDeviceToHost);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	// Free allocated memory on GPU
	if (cudaFree(d_tmp_mutant) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
    
	// Free allocated memory on GPU
	if (cudaFree(d_sequence) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	//printf("%s\n%s\n",sequence, mutant);
	
	return mutant;
}

__device__ void CUDAGetNK(int mutant_num, int seq2_len, int* n, int* k)
{
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

//	returns the score for a given mutant and offset, the func will compare 2 chars each loop, one from seq1 and another from seq2, it will skip the letters in indexes N and K,
//	that will make the comparison of the two strings a comparison between seq2 and a mutant sequence(N,K)
__device__ float calcMutantScore(char* seq1, char* seq2, float* d_conservative_matrix,int len2, int n, int k, int index, int offset)
{
	// N and K needs to be minus 1 because for thread i=0 we will get (N,K)=(1,2) and we don't want to skip 1,2 we want to skip indexes (0,1) therefore we subtract 1 from each
	n = n - 1;
	k = k - 1;
	float score = 0;
	int i = 0, j = i;
	for (i = 0; i < len2 - 2; i++, j++)
	{
		if(j == n)
			j++; // skipping the comparison with index n
		if(j == k)
			j++; // skipping the comparison with index k
		float tmp_score = d_conservative_matrix[(seq1[i] - 'A') * 26 + (seq2[j] - 'A')]; // the comparison matrix is used here as a 2d array but in fact it is a 1d array
		score += tmp_score;
	}	
	//printf("(%d,%d) %1.2f\n",n, k, score);
	return score;	
}

//	each thread is responsible for a single mutant and calculates its best score and offset, example:
//	thread i=0 will calculate the result for (N=1,K=2) and store the best score it could find in d_mutantsBestScores[0] and best offset it could find in d_mutantsBestOffsets[0]
__global__ void calcMutantBestScoreKernel(char* d_seq1, char* d_seq2, float* d_comp_matrix, float* d_mutantsBestScores, int* d_mutantsBestOffsets, int num_mutants, int maxOffset, int len2)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int n,k;
	float bestScore = -10000;
	int offset = 0;
	if (i < num_mutants)
	{
		CUDAGetNK(i+1, len2, &n, &k);
		//printf("i = %d, (%d,%d)\n", i, n, k);
		for (int j = 0; j < maxOffset; j++)
		{
			float score = calcMutantScore(&d_seq1[j], d_seq2, d_comp_matrix, len2, n, k, i, j);
			// asking if the score we found for the given (N,K) and offset is bigger than the best score we found so far
			if (score > bestScore)
			{
				bestScore = score;
				offset = j;	
			}
		}
		// storing the final results in the 2 vectors
		d_mutantsBestScores[i]	= bestScore;
		d_mutantsBestOffsets[i] = offset;
	}
}

/*
	this function calculates the best scores and offsets for all the mutants, the result of those results will be in two vectors: mutantsBestScores, mutantsBestOffsets. the
	length of those 2 vectors is the number of mutants for a given seq2.
	this function uses CUDA to calculate the score and offset for each mutant and stores the results in the previously mentioned vectors.
	later it uses OMP to get the best score, best offset and (N,K) for the given seq2.
	return the string result that is to be printed.
*/
char* calc_best_score_CUDA(char* seq1, char* seq2, float* comp_matrix){
	
	int seq1_len = strlen(seq1);
	int seq2_len = strlen(seq2);
	
	//printf("seq1_len = %d, seq2_len = %d\n", seq1_len, seq2_len);
	//printf("seq2 = %s\n", seq2);
	
	// calc maxmimum offset and number of mutants
	int maxOffset = seq1_len - (seq2_len - 2) + 1;
	int num_of_mutants = seq2_len * (seq2_len - 1) / 2;
	
	// allocate memory in CPU for 2 vectors
	float* mutantsBestScores = (float*) malloc(num_of_mutants * sizeof(float));
	int* mutantsBestOffsets = (int*) malloc(num_of_mutants * sizeof(int));
	
	// allocate d_seq1, d_seq2, d_comp_matrix memory and copy data to device
	char* d_seq1 = NULL;
	char* d_seq2 = NULL;
	float* d_comp_matrix = NULL;
	CUDA_MEM_INIT_COPY(d_seq1, seq1, seq1_len, char);
	CUDA_MEM_INIT_COPY(d_seq2, seq2, seq2_len, char);
	CUDA_MEM_INIT_COPY(d_comp_matrix, comp_matrix, SIZE_OF_COMP_MATRIX, float);
	
	// allocate memory for d_mutantsBestScores & d_mutantsBestOffsets and dont copy data! - the data will be filled in the device
	float* d_mutantsBestScores = NULL; 
	int* d_mutantsBestOffsets = NULL;
	CUDA_MEM_INIT(d_mutantsBestScores, num_of_mutants, float);
	CUDA_MEM_INIT(d_mutantsBestOffsets, num_of_mutants, int);
	
	int threads = 256;
	int blocks = (num_of_mutants + threads - 1) / threads;
	
	calcMutantBestScoreKernel<<<blocks, threads>>>(d_seq1, d_seq2, d_comp_matrix, d_mutantsBestScores, d_mutantsBestOffsets, num_of_mutants, maxOffset, seq2_len);
	
	cudaMemcpy(mutantsBestScores, d_mutantsBestScores, num_of_mutants * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(mutantsBestOffsets, d_mutantsBestOffsets, num_of_mutants * sizeof(int), cudaMemcpyDeviceToHost);
	for(int i = 0; i < num_of_mutants; i++){
		//printf("i = %d, score = %1.2f\n", i, mutantsBestScores[i]);
	}
	/*
	float maxScore = -10000;
	int bestOffset = 0;
	int bestMutantNum = -1;
	
	for (int i = 0; i < num_of_mutants; i++)
	{
		if (mutantsBestScores[i] > maxScore)
		{
			maxScore = mutantsBestScores[i];
			bestOffset = mutantsBestOffsets[i];
			bestMutantNum = i;
		}
	}
	int n,k;
	CPUGetNK(bestMutantNum + 1, seq2_len, &n, &k);
	printf("mutant num: %d, MS(%d,%d), score: %1.2f, offset: %d\n", bestMutantNum, n, k, maxScore, bestOffset);
	*/
	
	// using OMP to get the best results	
	char* final_result = calcBestScoreOmp(mutantsBestScores, mutantsBestOffsets, num_of_mutants, seq2_len);
	printf("%s", final_result);
	return final_result;
}
