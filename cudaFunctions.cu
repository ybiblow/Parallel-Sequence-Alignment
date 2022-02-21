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

__global__ void incrementByOne(int *arr, int numElements) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // Increment the proper value of the arrray according to thread ID 
    if (i < numElements)
        arr[i]++;
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
			
    /*
	//printf("getting mutant (%d,%d) length = %d\n", m, n, len);	
	int i;
	int f_index = m - 1;
	int e_index = n - 1;
	char* mutant = (char*)malloc((len - 1) * sizeof(char));	
	#pragma omp parallel for shared(f_index, e_index, mutant)
	for(i = 0; i < len; i++){
		if(i < f_index)
			mutant[i] = sequence[i];
		else if(i > f_index && i < e_index)
			mutant[i - 1] = sequence[i];
		else if(i > e_index)
			mutant[i - 2] = sequence[i];
	}
	mutant[len - 2] = '\0';	
	return mutant;
	*/


    // Increment the proper value of the arrray according to thread ID 
}

int computeOnGPU(int *data, int numElements) {
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    size_t size = numElements * sizeof(float);
  

    // Allocate memory on GPU to copy the data from the host
    int *d_A;
    err = cudaMalloc((void **)&d_A, size);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy data from host to the GPU memory
    err = cudaMemcpy(d_A, data, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }


    // Launch the Kernel
    int threadsPerBlock = 256;
    int blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
    incrementByOne<<<blocksPerGrid, threadsPerBlock>>>(d_A, numElements);
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to launch vectorAdd kernel -  %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the  result from GPU to the host memory.
    err = cudaMemcpy(data, d_A, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Free allocated memory on GPU
    if (cudaFree(d_A) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    return 0;
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

void calc_best_score_CUDA(char* seq1, char* mutant, float* weights,int* best_offset, float* best_score){
	int i;
	int offset = strlen(seq1) - strlen(mutant);
	int tmp_score;
	*best_score = 0;
	//printf("offset = %d\n", offset);
	for(i = 0; i <= offset; i++){
		// working on changing calc_similarity_CUDA to work with cuda
		char* similarity = calc_similarity_CUDA(seq1, mutant, i);
		tmp_score = calc_score(similarity, strlen(mutant), &weights[0]);
		if(tmp_score > *best_score){
			*best_score = tmp_score;
			*best_offset = i;
		}
		free(similarity);
	}
	//printf("Best offset is: %d\n", *best_offset);
	
}

char* calc_similarity_CUDA(char* seq1, char* mutant, int offset){
	int len = strlen(mutant);
	size_t size = len * sizeof(char);
	char* result = (char*)malloc(len * sizeof(char));
	
	// Error code to check return values for CUDA calls
	cudaError_t err = cudaSuccess;	
  
	// Allocate memory on GPU to copy the data from the host
	char* d_seq1;
	err = cudaMalloc((void **)&d_seq1, strlen(seq1) * sizeof(char));      
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Allocate memory on GPU to copy the data from the host
	char* d_mutant;
	err = cudaMalloc((void **)&d_mutant, size);      
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Allocate memory on GPU to copy the data from the host
	char* d_result;
	err = cudaMalloc((void **)&d_result, size);      
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Copy data from host to the GPU memory
	err = cudaMemcpy(d_seq1, seq1, strlen(seq1) * sizeof(char), cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Copy data from host to the GPU memory
	err = cudaMemcpy(d_mutant, mutant, size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Launch the Kernel
	int threadsPerBlock = 256;
	int blocksPerGrid =(len + threadsPerBlock - 1) / threadsPerBlock;
	calc_similarity_CUDA_Kernel<<<blocksPerGrid, threadsPerBlock>>>(d_seq1, d_mutant, len, d_result, offset);
	err = cudaGetLastError();
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to launch vectorAdd kernel -  %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Copy the  result from GPU to the host memory.
	err = cudaMemcpy(result, d_result, size, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Free allocated memory on GPU
	if (cudaFree(d_seq1) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Free allocated memory on GPU
	if (cudaFree(d_mutant) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	// Free allocated memory on GPU
	if (cudaFree(d_result) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	printf("offset = %d\n", offset);
	return result;
}
