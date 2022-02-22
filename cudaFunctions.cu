#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "myProto.h"

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

void calc_best_score_CUDA(char* seq1, char* mutant, float* weights, int* best_offset, float* best_score){
	int seq1_len = strlen(seq1);
	size_t seq1_size = seq1_len * sizeof(char);
	int mutant_len = strlen(mutant);
	size_t mutant_size = mutant_len * sizeof(char);
	int offset = seq1_len - mutant_len;
	size_t d_score_size = (offset + 1) * sizeof(float);
	size_t d_weights_size = 4 * sizeof(float);
	
	float* score_arr = (float*)malloc((offset + 1) * sizeof(float));
	// Error code to check return values for CUDA calls
	cudaError_t err = cudaSuccess;
	
	// Allocate memory on GPU to copy the data from the host - d_seq1
	char* d_seq1;
	err = cudaMalloc((void **)&d_seq1, seq1_size);      
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Allocate memory on GPU to copy the data from the host - d_mutant
	char* d_mutant;
	err = cudaMalloc((void **)&d_mutant, mutant_size);      
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Allocate memory on GPU to copy the data from the host - d_score_arr
	float* d_score_arr;
	err = cudaMalloc((void **)&d_score_arr, d_score_size);      
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Allocate memory on GPU to copy the data from the host - d_weights
	float* d_weights;
	err = cudaMalloc((void **)&d_weights, d_weights_size);      
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Copy data from host to the GPU memory - d_seq1
	err = cudaMemcpy(d_seq1, seq1, seq1_size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Copy data from host to the GPU memory - d_mutant
	err = cudaMemcpy(d_mutant, mutant, mutant_size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Copy data from host to the GPU memory - d_weights
	err = cudaMemcpy(d_weights, weights, d_weights_size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	//void calcBestScore(char* d_seq1, char* d_mutant, int* d_score_arr, float* weights, int offset, int mutant_len)	

	// Launch the Kernel
	int threadsPerBlock = 256;
	int blocksPerGrid =(offset + threadsPerBlock - 1) / threadsPerBlock;
	calcBestScore<<<blocksPerGrid, threadsPerBlock>>>(d_seq1, d_mutant, d_score_arr, d_weights, offset, mutant_len);
	err = cudaGetLastError();
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to launch vectorAdd kernel -  %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Copy the  result from GPU to the host memory. d_score_arr ---> score_arr
	err = cudaMemcpy(score_arr, d_score_arr, d_score_size, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// score_arr has the info we need to get the best score
	int score = -1;
	*best_offset = -1;
	int i;
	for(i = 0; i <= offset; i++){
		if(score_arr[i] > score){
			score = score_arr[i];
			*best_offset = i;
		}
	}
	*best_offset = offset;
	*best_score = score;
	
	// Free allocated memory on GPU - d_seq1
	if (cudaFree(d_seq1) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Free allocated memory on GPU - d_mutant
	if (cudaFree(d_mutant) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	// Free allocated memory on GPU - d_score_arr
	if (cudaFree(d_score_arr) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Free allocated memory on GPU - d_weights
	if (cudaFree(d_weights) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	free(score_arr);
	//*best_score = ;
	//*best_offset = ;
	//printf("I'm here, This is where you want to be!!!\n");
}
