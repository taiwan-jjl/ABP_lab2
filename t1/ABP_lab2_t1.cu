
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdbool.h>
#include <time.h>


__device__ float reduction(float *sdata
                           )
{
    float sum = 0;
    atomicAdd(&sum, sdata);
    __syncthreads();
    return sum;
}




__global__ void compute(uint   M,
                        uint   N,
                        float *x,
                        float *y)
{
    const int idx = threadIdx.x + blockIdx.x * blockDim.x;
    extern __shared__ float sdata[];

    if(idx < N){
        sdata[threadIdx.x] = x[idx] * y[idx];
    }
}




int main(int argc, char *argv[]){

    // hardware parameter for A3000
    //const int num_SM = 32;
    //const int num_warp = 32;

    // input parameter
    uint M = 1024;
    uint N = 512;
    size_t num_thread = 512;
    dim3 num_block = ( (uint)(ceilf((float)N / (float)num_thread) ), M );
    
    float *A, *x, *y;
    // allocate memory on the device
    cudaMalloc( &A, M*N * sizeof(float));
    cudaMalloc( &x, N * sizeof(float));
    cudaMalloc( &y, N * sizeof(float));

    // initialize variable value
    cudaMemset( A, 1, M*N);
    cudaMemset( x, 3, N);
    cudaMemset( y, 5, N);
    cudaDeviceSynchronize();

    // compute
    compute<<<num_block, num_thread, num_thread*sizeof(float), 0>>>(M, N, x, y );


    // Free the memory on the device
    cudaFree(A);
    cudaFree(x);
    cudaFree(y);

    return 0;
}