#define HOST_EMULATION
#include "cuda_runtime_api.h"
float atomicAdd(float *a, float b) { float o = *a; *a += b; return o; }
double atomicAdd(double *a, double b) { double o = *a; *a += b; return o; }
int atomicAdd(int *a, int b) { int o = *a; *a += b; return o; }
unsigned atomicAdd(unsigned *a, unsigned b) { unsigned o = *a; *a += b; return o; }
unsigned long long atomicAdd(unsigned long long *a, unsigned long long b) { auto o = *a; *a += b; return o; }
int atomicExch(int *a, int b) { int o = *a; *a = b; return o; }
int atomicCAS(int *a, int c, int v) { int o = *a; if(o == c) *a = v; return o; }
unsigned long long atomicCAS(unsigned long long *a, unsigned long long c, unsigned long long v) { auto o = *a; if(o == c) *a = v; return o; }
uint3 threadIdx, blockIdx; dim3 blockDim, gridDim;
