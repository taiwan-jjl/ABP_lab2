
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdbool.h>
#include <time.h>


__global__ void compute(const int   M,
                        const int   N,
                        const float *x,
                        const float *y,
                        const float *z)
{
    const int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if(idx < N){
        z[idx] = a * x[idx] + y[idx];
    }
}




int main(int argc, char *argv[]){

    // hardware parameter for A3000
    const int num_SM = 32;
    const int num_warp = 32;

    // input parameter
    int M = 1024;
    int N = 512;
    int num_thread = 1024;
    int num_block = (int)(ceilf((float)M*N / (float)num_thread))
    
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
    compute<<<num_block, num_thread>>>()


    // Free the memory on the device
    cudaFree(A);
    cudaFree(x);
    cudaFree(y);

    return 0;
}