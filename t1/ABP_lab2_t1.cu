
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdbool.h>
#include <time.h>

#define num_thread 512


__global__ void result(float *y,
                       uint count
                           )
{
    for(int i=0;i<count;i++){
        printf("i=%i, %f\n",i, y[i]);
    }
    
}




__global__ void compute(uint   M,
                        uint   N,
                        float *A,
                        float *x,
                        float *y)
{
    extern __shared__ float sdata[];
    const int idx = blockIdx.y * (int)N + blockIdx.x * blockDim.x + threadIdx.x;
    const int idxx = threadIdx.x + blockIdx.x * blockDim.x;

    if(idxx < N){
        sdata[threadIdx.x] = A[idx] * x[idxx];
    }
    __syncthreads();

    if(threadIdx.x != 0){
        atomicAdd(&sdata[0], sdata[threadIdx.x]);
    }
    __syncthreads();
    if(threadIdx.x == 0){
        atomicAdd(&y[blockIdx.y], sdata[0]);
    }
}




__global__ void setmemoryf(float *A,
                           float value, 
                           size_t count)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx<count){
        A[idx] = value;
    }
    // if(threadIdx.x==0 && blockIdx.x==0){
    // for(int i=0;i<count;i++){
    //     printf("i=%i, %f\n",i, A[i]);
    // }}
}




int main(int argc, char *argv[]){

    // hardware parameter for A3000
    //const int num_SM = 32;
    //const int num_warp = 32;
    // 4 warp schedulers.

    // input parameter
    uint M = 1024;
    uint N = 1024;
    //size_t num_thread = 512;
    dim3 grid( (uint)(ceilf((float)N / (float)num_thread) ), M );
    
    float *A, *x, *y;
    // allocate memory on the device
    cudaMalloc( &A, M*N * sizeof(float));
    cudaMalloc( &x, N * sizeof(float));
    cudaMalloc( &y, M * sizeof(float));

    // initialize variable value
    //cudaMemset is only for integer!!!
    //cudaMemset( A, 1, M*N);
    //cudaMemset( x, 3, N);
    //cudaMemset( y, 5, M);
    setmemoryf<<<(uint)(ceilf((float)M*N / (float)num_thread) ), num_thread, 0, 0>>>(A, 1.0f, M*N);
    setmemoryf<<<(uint)(ceilf((float)N / (float)num_thread) ), num_thread, 0, 0>>>(x, 3.0f, N);
    setmemoryf<<<(uint)(ceilf((float)M / (float)num_thread) ), num_thread, 0, 0>>>(y, 0.0f, M);
    cudaDeviceSynchronize();

    // compute and timing
    struct timespec start, end;
    timespec_get(&start, TIME_UTC);

    int repeat = 100000;
    for(int i=0;i<repeat;i++){
        compute<<<grid, num_thread, num_thread*sizeof(float), 0>>>(M, N, A, x, y );
        cudaDeviceSynchronize();
    }

    timespec_get(&end, TIME_UTC);
    time_t d_sec  = end.tv_sec  - start.tv_sec;
    long   d_nsec = end.tv_nsec - start.tv_nsec;
    double total_time = (double)d_sec + (double)d_nsec/1000000000.0;
    printf("time=%f\n", total_time);
    double perf = 

    // show result
    //result<<<1,1>>>(y,M);


    // Free the memory on the device
    cudaFree(A);
    cudaFree(x);
    cudaFree(y);

    return 0;
}