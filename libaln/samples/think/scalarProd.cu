/* This code is partly from NVIDIA CUDA TOOLKIT.
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

/*
 * This sample calculates scalar products of a
 * given set of input vector pairs
 */



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <helper_functions.h>
#include <helper_cuda.h>


void callCUDA();
void initGPU();
void execGPU();
void CheckResult();

///////////////////////////////////////////////////////////////////////////////
// Calculate scalar products of VectorN vectors of ElementN elements on CPU
///////////////////////////////////////////////////////////////////////////////
extern "C"
void scalarProdCPU(
    float *h_C,
    float *h_A,
    float *h_B,
    int vectorN,
    int elementN
);



///////////////////////////////////////////////////////////////////////////////
// Calculate scalar products of VectorN vectors of ElementN elements on GPU
///////////////////////////////////////////////////////////////////////////////
#include "scalarProd_kernel.cuh"



////////////////////////////////////////////////////////////////////////////////
// Helper function, returning uniformly distributed
// random float in [low, high] range
////////////////////////////////////////////////////////////////////////////////
//float RandFloat(float low, float high)
//{
//    float t = (float)rand() / (float)RAND_MAX;
//    return (1.0f - t) * low + t * high;
//}



///////////////////////////////////////////////////////////////////////////////
// Data configuration
///////////////////////////////////////////////////////////////////////////////

//Total number of input vector pairs; arbitrary
//Number of elements per vector; arbitrary,
//but strongly preferred to be a multiple of warp size
//to meet memory coalescing constraints

extern int ELEMENT_N;
extern int nNumberLFNs;
#define VECTOR_N nNumberLFNs
//Total number of data elements
const int    DATA_N = VECTOR_N * ELEMENT_N;

const int   DATA_SZ = DATA_N * sizeof(float);
const int RESULT_SZ = VECTOR_N  * sizeof(float);
float* h_C_CPU, * h_C_GPU;
float* d_A, * d_B, * d_C;

void callCUDA()
{
}

///////////////////////////////////////////////////////////////////////////////
// Main program
///////////////////////////////////////////////////////////////////////////////
void initGPU()
{
    int i;

    printf("%s Starting...\n\n", argv[0]);

    // use command-line specified CUDA device, otherwise use device with highest Gflops/s
    findCudaDevice(argc, (const char**)argv);

   

    printf("Initializing data...\n");
    printf("...allocating CPU memory.\n");
    // h_A = (float*)malloc(DATA_SZ); // these two are allocated in main, and are of maximal size
    //h_B = (float*)malloc(DATA_SZ);
    h_C_CPU = (float*)malloc(RESULT_SZ);
    h_C_GPU = (float*)malloc(RESULT_SZ);

    printf("...allocating GPU memory.\n");
    checkCudaErrors(cudaMalloc((void**)&d_A, DATA_SZ));
    checkCudaErrors(cudaMalloc((void**)&d_B, DATA_SZ));
    checkCudaErrors(cudaMalloc((void**)&d_C, RESULT_SZ));

    execGPU();
    // We check inside a subroutine
    // CheckResult(h_C_CPU, h_C_GPU, h_A, h_B, VECTOR_N, ELEMENT_N);
    CheckResult();
    free(h_C_GPU);
    free(h_C_CPU);
    free(h_B);
    free(h_A);


    printf("Shutting down...\n");
    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_B));
    checkCudaErrors(cudaFree(d_A));

}

void execGPU()
{
    StopWatchInterface* hTimer = NULL;
    sdkCreateTimer(&hTimer);
    printf("...copying input data to GPU mem.\n");
    //Copy options data to GPU memory for further processing
    checkCudaErrors(cudaMemcpy(d_A, h_A, DATA_SZ, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_B, h_B, DATA_SZ, cudaMemcpyHostToDevice));
    printf("Data init done.\n");


    printf("Executing GPU kernel...\n");
    checkCudaErrors(cudaDeviceSynchronize());
    sdkResetTimer(&hTimer);
    sdkStartTimer(&hTimer);
    scalarProdGPU << <128, 256 >> > (d_C, d_A, d_B, VECTOR_N, ELEMENT_N);
    getLastCudaError("scalarProdGPU() execution failed\n");
    checkCudaErrors(cudaDeviceSynchronize());
    sdkStopTimer(&hTimer);
    printf("GPU time: %f msecs.\n", sdkGetTimerValue(&hTimer));
    sdkDeleteTimer(&hTimer);

    printf("Reading back GPU result...\n");
    //Read back GPU results to compare them to CPU results
    checkCudaErrors(cudaMemcpy(h_C_GPU, d_C, RESULT_SZ, cudaMemcpyDeviceToHost));
}


void CheckResult()
{
    double delta, ref, sum_delta, sum_ref, L1norm;

    printf("Checking GPU results...\n");
    printf("..running CPU scalar product calculation\n");
    scalarProdCPU(h_C_CPU, h_A, h_B, VECTOR_N, ELEMENT_N);

    printf("...comparing the results\n");
    //Calculate max absolute difference and L1 distance
    //between CPU and GPU results
    sum_delta = 0;
    sum_ref   = 0;

    for (int i = 0; i < VECTOR_N; i++)
    {
        delta = fabs(h_C_GPU[i] - h_C_CPU[i]);
        ref   = h_C_CPU[i];
        sum_delta += delta;
        sum_ref   += ref;
    }

    L1norm = sum_delta / sum_ref;



    printf("L1 error: %E\n", L1norm);
    printf((L1norm < 1e-6) ? "Test passed\n" : "Test failed!\n");
    exit(L1norm < 1e-6 ? EXIT_SUCCESS : EXIT_FAILURE);
}
