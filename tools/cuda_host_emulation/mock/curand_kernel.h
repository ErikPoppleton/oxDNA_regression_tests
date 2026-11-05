#pragma once
#include "cuda_runtime_api.h"
struct curandStateXORWOW { int d; }; typedef curandStateXORWOW curandState;
__device__ float curand_normal(curandState*); __device__ double curand_normal_double(curandState*); __device__ float curand_uniform(curandState*); __device__ void curand_init(unsigned long long, unsigned long long, unsigned long long, curandState*);
