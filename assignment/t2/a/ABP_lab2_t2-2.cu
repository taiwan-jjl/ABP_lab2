
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdbool.h>
#include <time.h>
#include "cublas_v2.h"

#define num_thread 512
#define warp 32


__global__ void result(float *y,
                       int count
                           )
{
    for(int i=0;i<count;i++){
        printf("i=%i, %f\n",i, y[i]);
    }
    
}




__global__ void compute(int    M,
                        int    N,
                        float *A,
                        float *x,
                        float *y)
{
    __shared__ float sx[warp];
    __shared__ float sy[warp];
    if(threadIdx.y==0){
        sy[threadIdx.x]=0.0f;
    }
               float temp=0.0f;

    const int mat_x = blockIdx.y*blockDim.y + threadIdx.y;
    const int mat_y = blockIdx.x*blockDim.x + threadIdx.x;
    const int mat_idx = mat_x*M + mat_y;

    if(threadIdx.x==0 && mat_x<N){
        sx[threadIdx.y] = x[mat_x];
    }
    __syncthreads();

    if(mat_x<N && mat_y<M){
        temp = A[mat_idx] * sx[threadIdx.y];
        atomicAdd(&sy[threadIdx.x], temp);
    }
    __syncthreads();

    if(threadIdx.y==0 && mat_x<N && mat_y<M){
        atomicAdd(&y[mat_y], sy[threadIdx.x]);
    }
    __syncthreads();

    // if(threadIdx.x==0 && threadIdx.y==0){
    //     printf("blk-x=%u, blk-y=%u \n", blockIdx.x, blockIdx.y);
    // }
}




__global__ void setmemoryf(float *A,
                           float value, 
                           int   count)
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
    // const int num_SM = 32;
    // const int num_warp = 32;
    // 4 warp schedulers.

    // input parameter
for(int test = 128; test<10000; test=test+32){

    int M = 10240;
    int N = 10240;
M=test;
N=test;

    //size_t num_thread = 512;
    //dim3 threadblock(warp, warp);
    //dim3 blockgrid( (unsigned int)(ceilf((float)M / (float)warp) ), (unsigned int)(ceilf((float)N / (float)warp) ) );
    
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
    setmemoryf<<<(unsigned int)(ceilf((float)M*N / (float)num_thread) ), num_thread, 0, 0>>>(A, 1.0f, M*N);
    setmemoryf<<<(unsigned int)(ceilf((float)N / (float)num_thread) ), num_thread, 0, 0>>>(x, 3.0f, N);
    setmemoryf<<<(unsigned int)(ceilf((float)M / (float)num_thread) ), num_thread, 0, 0>>>(y, 0.0f, M);
    cudaDeviceSynchronize();

    // compute and timing

        cublasHandle_t handle;
        cublasStatus_t stat = cublasCreate(&handle);
        if (stat != CUBLAS_STATUS_SUCCESS){
            printf("CUBLAS initialization failed\n");
            return 1;
        }

        const float alpha = 1.0f;
        const float beta  = 0.0f;




    int repeat = 1000;

    struct timespec start, end;
    timespec_get(&start, TIME_UTC);


    for(int i=0;i<repeat;i++){


        stat =cublasSgemv(handle, CUBLAS_OP_N, 
                          M, N, 
                          &alpha, A, M, 
                          x, 1, 
                          &beta, 
                          y, 1);
        cudaDeviceSynchronize();


    }

    timespec_get(&end, TIME_UTC);

        if(stat != CUBLAS_STATUS_SUCCESS){
            printf("CUBLAS operation failed\n");
            return 1;
        }
        
        cublasDestroy(handle);






    time_t d_sec  = end.tv_sec  - start.tv_sec;
    long   d_nsec = end.tv_nsec - start.tv_nsec;
    double total_time = (double)d_sec + (double)d_nsec/1000000000.0;
    //printf("time=%f\n", total_time);
    double perf = 1.0e-9 * ((double)M*(double)N + (double)M + (double)N) * 4.0 / total_time * (double)repeat;
    printf("test= %i , time(sec)= %f ,memory throughput = %f GByte/s\n", test, total_time, perf);

    // show result
    //result<<<1,1>>>(y,M);


    // Free the memory on the device
    cudaFree(A);
    cudaFree(x);
    cudaFree(y);
}
    return 0;
}