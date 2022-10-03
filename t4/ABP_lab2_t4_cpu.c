
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdbool.h>
#include <time.h>

#define num_thread 512
#define warp 32


void result(float *y,
                       int count
                           )
{
    for(int i=0;i<count;i++){
        printf("i=%i, %f\n",i, y[i]);
    }
    
}




void setmemoryf(float *A,
                           float value, 
                           int   count)
{
    for(int i=0; i<count; i++){
        A[i] = value;
    }
}




int main(int argc, char *argv[]){

    // hardware parameter for A3000
    // const int num_SM = 32;
    // const int num_warp = 32;
    // 4 warp schedulers.

    // input parameter
//for(int test = 128; test<10000; test=test+32){

    int M = 3;
    int N = 3;
    int K = 3;
//M=test;
//N=test;
//K=test;
    
    // allocate memory on the device
    float *A = malloc(M*K * sizeof(float));
    float *x = malloc(K*N * sizeof(float));
    float *y = malloc(M*N * sizeof(float));

    // initialize variable value
    setmemoryf(A, 1.0f, M*K);
    setmemoryf(x, 3.0f, K*N);
    setmemoryf(y, 0.0f, M*N);

    // compute and timing
    struct timespec start, end;
    timespec_get(&start, TIME_UTC);

    int repeat = 1;
    for(int i=0;i<repeat;i++){

        for(int idxN=0; idxN<N; idxN++){
            for(int idxM=0; idxM<M; idxM++){
                for(int idxK=0; idxK<K; K++){

                    y[idxN*M + idxM] = y[idxN*M + idxM] + A[idxM + idxK*M] * x[idxN*M + idxK];

                }
            }
        }

    }

    timespec_get(&end, TIME_UTC);
    time_t d_sec  = end.tv_sec  - start.tv_sec;
    long   d_nsec = end.tv_nsec - start.tv_nsec;
    double total_time = (double)d_sec + (double)d_nsec/1000000000.0;
    //printf("time=%f\n", total_time);
    double perf = 1.0e-9 * ((double)M*(double)N*(double)K*(double)2) / total_time * (double)repeat;
    //printf("test= %i , time(sec)= %f ,memory throughput = %f GFlop/s\n", test, total_time, perf);

    // show result
    result(y,M*N);


    // Free the memory on the device
    free(A);
    free(x);
    free(y);
//}
    return 0;
}