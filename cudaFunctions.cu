#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "myProto.h"

__global__ void incrementByOne(int *arr, int numElements) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // Increment the proper value of the arrray according to thread ID 
    if (i < numElements)
        arr[i]++;
}

__global__ void get_Mutant_CUDA_Kernel(char *sequence, char* tmp_mutant, int len, int m, int n) {
		int i = blockDim.x * blockIdx.x + threadIdx.x;
		if (i < len){
			if(i < (m - 1))
				tmp_mutant[i] = sequence[i];
			else if(i > (m - 1) && i < (n - 1))
				tmp_mutant[i - 1] = sequence[i];
			else if(i > (n - 1))
				tmp_mutant[i - 2] = sequence[i];
		}
		if(i == len - 1)
			tmp_mutant[i - 2] = '\0';
			
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
	
	// Error code to check return values for CUDA calls
	cudaError_t err = cudaSuccess;

	size_t size = len * sizeof(char);
	char* mutant = (char*) malloc(size);

	// Allocate memory on GPU to copy the data from the host
	char* *d_sequence;
	err = cudaMalloc((void **)&d_sequence, size);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	size_t size_tmp_mutant = (len - 2) * sizeof(char);
	char* *d_tmp_mutant;

	// Allocate memory on GPU to copy the data from the host
	err = cudaMalloc((void **)&d_tmp_mutant, size_tmp_mutant);
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
	err = cudaMemcpy(mutant, d_sequence, size, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	
	
	// Free allocated memory on GPU
	if (cudaFree(d_sequence) != cudaSuccess) {
		fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	//mutant[len] = '\0';
	return mutant;
}

