#pragma once
#include <cstddef>
#include <cmath>
#ifdef HOST_EMULATION
#define __restrict__
#define __host__
#define __device__
#define __global__
#define __constant__
#define __shared__
#define __forceinline__ inline
#define __align__(n) __attribute__((aligned(n)))
#define __builtin_align__(n) __attribute__((aligned(n)))
#else
#define __restrict__
#define __host__ __attribute__((host))
#define __device__ __attribute__((device))
#define __global__ __attribute__((global))
#define __constant__ __attribute__((constant))
#define __shared__ __attribute__((shared))
#define __forceinline__ __inline__ __attribute__((always_inline))
#define __align__(n) __attribute__((aligned(n)))
#define __builtin_align__(n) __attribute__((aligned(n)))
#endif
#define CUDART_VERSION 12000
typedef int cudaError_t; enum { cudaSuccess = 0 };
typedef int cudaStream_t; typedef int CUresult; typedef int cudaEvent_t;
enum cudaMemcpyKind { cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost, cudaMemcpyDeviceToDevice, cudaMemcpyHostToHost };
struct dim3 { unsigned x, y, z; __host__ __device__ dim3(unsigned a = 1, unsigned b = 1, unsigned c = 1) : x(a), y(b), z(c) {} };
struct uint3 { unsigned x, y, z; };
extern __device__ uint3 threadIdx, blockIdx; extern __device__ dim3 blockDim, gridDim;
struct __align__(16) float4 { float x, y, z, w; };
struct float3 { float x, y, z; };
struct float2 { float x, y; };
struct __align__(16) double4 { double x, y, z, w; };
struct double3 { double x, y, z; };
struct int4 { int x, y, z, w; }; struct int2 { int x, y; };
struct uint4 { unsigned x,y,z,w; };
struct uchar4 { unsigned char x,y,z,w; };
__host__ __device__ inline float4 make_float4(float a, float b, float c, float d) { return {a, b, c, d}; }
__host__ __device__ inline float3 make_float3(float a, float b, float c) { return {a, b, c}; }
__host__ __device__ inline double4 make_double4(double a, double b, double c, double d) { return {a, b, c, d}; }
__host__ __device__ inline int4 make_int4(int a, int b, int c, int d) { return {a, b, c, d}; }
__device__ float atomicAdd(float *, float); __device__ double atomicAdd(double *, double); __device__ int atomicAdd(int *, int); __device__ unsigned atomicAdd(unsigned *, unsigned);
__device__ unsigned long long atomicAdd(unsigned long long *, unsigned long long);
__device__ int atomicExch(int*, int); __device__ int atomicCAS(int*,int,int); __device__ unsigned long long atomicCAS(unsigned long long*, unsigned long long, unsigned long long);
__device__ long long __double_as_longlong(double); __device__ double __longlong_as_double(long long);
__device__ void __syncthreads();
#ifndef HOST_EMULATION
__device__ float __int_as_float(int); __device__ int __float_as_int(float);
#endif
#ifndef HOST_EMULATION
template<class T> cudaError_t cudaMemcpyToSymbol(const T &, const void *, size_t, size_t = 0, cudaMemcpyKind = cudaMemcpyHostToDevice);
#endif
template<class T> cudaError_t cudaMemcpyFromSymbol(void *, const T &, size_t, size_t = 0, cudaMemcpyKind = cudaMemcpyDeviceToHost);
#ifndef HOST_EMULATION
cudaError_t cudaMemcpy(void *, const void *, size_t, cudaMemcpyKind);
#endif
 cudaError_t cudaMemset(void *, int, size_t);
#ifndef HOST_EMULATION
cudaError_t cudaMalloc(void **, size_t);
#endif
#ifndef HOST_EMULATION
template<class T> cudaError_t cudaMalloc(T **p, size_t s) { return cudaMalloc((void**)p, s); }
#endif
cudaError_t cudaMallocManaged(void **, size_t, unsigned = 1); template<class T> cudaError_t cudaMallocManaged(T **p, size_t s, unsigned f = 1) { return cudaMallocManaged((void**)p, s, f); }
#ifndef HOST_EMULATION
cudaError_t cudaFree(void *);
#endif
cudaError_t cudaGetLastError(); const char *cudaGetErrorString(cudaError_t); cudaError_t cudaDeviceSynchronize(); cudaError_t cudaThreadSynchronize();
cudaError_t cudaMallocHost(void **, size_t); cudaError_t cudaFreeHost(void*);
struct cudaDeviceProp { char name[256]; size_t totalGlobalMem; int major, minor, multiProcessorCount, maxThreadsPerBlock; int maxGridSize[3]; int maxThreadsDim[3]; int warpSize; size_t sharedMemPerBlock; int regsPerBlock; int clockRate; int memoryClockRate; int memoryBusWidth;};
cudaError_t cudaGetDeviceProperties(cudaDeviceProp *, int); cudaError_t cudaGetDeviceCount(int*); cudaError_t cudaSetDevice(int); cudaError_t cudaGetDevice(int*);
extern "C" int __cudaPushCallConfiguration(dim3, dim3, size_t = 0, void * = 0);
cudaError_t cudaLaunchKernel(const void*, dim3, dim3, void**, size_t, cudaStream_t);
#define CUDA_SAFE_CALL_NO_SYNC(x) x
typedef unsigned long long cudaTextureObject_t;
struct cudaChannelFormatDesc { int x, y, z, w; int f; };
template<class T> cudaChannelFormatDesc cudaCreateChannelDesc();
enum cudaResourceType { cudaResourceTypeLinear };
struct cudaResourceDesc { cudaResourceType resType; union { struct { void *devPtr; cudaChannelFormatDesc desc; size_t sizeInBytes; } linear; } res; };
enum cudaTextureReadMode { cudaReadModeElementType };
struct cudaTextureDesc { cudaTextureReadMode readMode; int addressMode[3]; int filterMode; int normalizedCoords; };
cudaError_t cudaCreateTextureObject(cudaTextureObject_t *, const cudaResourceDesc *, const cudaTextureDesc *, const void *);
cudaError_t cudaDestroyTextureObject(cudaTextureObject_t);
template<class T> __device__ T tex1Dfetch(cudaTextureObject_t, int);
struct int3 { int x, y, z; };
__host__ __device__ inline int3 make_int3(int a, int b, int c) { return {a, b, c}; }
#ifndef HOST_EMULATION
__device__ float sqrtf(float); __device__ double sqrt(double); __device__ float sqrt(float);
__device__ float rintf(float); __device__ double rint(double); __device__ float rint(float);
__device__ float floorf(float); __device__ double floor(double); __device__ float floor(float);
__device__ float ceilf(float); __device__ double ceil(double); __device__ float ceil(float);
__device__ float fabsf(float); __device__ double fabs(double); __device__ float fabs(float);
__device__ float expf(float); __device__ double exp(double); __device__ float exp(float);
__device__ float logf(float); __device__ double log(double); __device__ float log(float);
__device__ float sinf(float); __device__ double sin(double); __device__ float sin(float);
__device__ float cosf(float); __device__ double cos(double); __device__ float cos(float);
__device__ float tanf(float); __device__ double tan(double); __device__ float tan(float);
__device__ float acosf(float); __device__ double acos(double); __device__ float acos(float);
__device__ float asinf(float); __device__ double asin(double); __device__ float asin(float);
__device__ float atanf(float); __device__ double atan(double); __device__ float atan(float);
__device__ float tanhf(float); __device__ double tanh(double); __device__ float tanh(float);
__device__ float coshf(float); __device__ double cosh(double); __device__ float cosh(float);
__device__ float sinhf(float); __device__ double sinh(double); __device__ float sinh(float);
__device__ float erfcf(float); __device__ double erfc(double); __device__ float erfc(float);
__device__ float erfcxf(float); __device__ double erfcx(double); __device__ float erfcx(float);
__device__ float erff(float); __device__ double erf(double); __device__ float erf(float);
__device__ float cbrtf(float); __device__ double cbrt(double); __device__ float cbrt(float);
__device__ float rsqrtf(float); __device__ double rsqrt(double); __device__ float rsqrt(float);
__device__ float log1pf(float); __device__ double log1p(double); __device__ float log1p(float);
__device__ float expm1f(float); __device__ double expm1(double); __device__ float expm1(float);
__device__ float roundf(float); __device__ double round(double); __device__ float round(float);
__device__ float truncf(float); __device__ double trunc(double); __device__ float trunc(float);
__device__ float exp2f(float); __device__ double exp2(double); __device__ float exp2(float);
__device__ float log2f(float); __device__ double log2(double); __device__ float log2(float);
__device__ float log10f(float); __device__ double log10(double); __device__ float log10(float);
__device__ float powf(float, float); __device__ double pow(double, double); __device__ float pow(float, float);
__device__ float atan2f(float, float); __device__ double atan2(double, double); __device__ float atan2(float, float);
__device__ float copysignf(float, float); __device__ double copysign(double, double); __device__ float copysign(float, float);
__device__ float fminf(float, float); __device__ double fmin(double, double); __device__ float fmin(float, float);
__device__ float fmaxf(float, float); __device__ double fmax(double, double); __device__ float fmax(float, float);
__device__ float fmodf(float, float); __device__ double fmod(double, double); __device__ float fmod(float, float);
__device__ float hypotf(float, float); __device__ double hypot(double, double); __device__ float hypot(float, float);
__device__ void sincosf(float, float*, float*); __device__ void sincos(double, double*, double*);
__device__ int min(int,int); __device__ int max(int,int); __device__ float min(float,float); __device__ float max(float,float); __device__ double min(double,double); __device__ double max(double,double);
__device__ int abs(int);
__device__ float __fdividef(float,float); __device__ float __expf(float); __device__ float __powf(float,float); __device__ float __logf(float); __device__ float __sinf(float); __device__ float __cosf(float);
__device__ int __clz(int); __device__ int __ffs(int); __device__ int __popc(unsigned);
__device__ unsigned __ballot_sync(unsigned, int); __device__ int __shfl_sync(unsigned,int,int); __device__ float __shfl_down_sync(unsigned, float, int);
__device__ void __threadfence();
#endif
cudaError_t cudaConfigureCall(dim3 gridSize, dim3 blockSize, size_t sharedSize = 0, cudaStream_t stream = 0);
inline void getLastCudaError(const char *) {}
template<class T> __device__ T __ldg(const T*);
#ifdef HOST_EMULATION
inline float __cosf(float x) { return cosf(x); }
inline float __sinf(float x) { return sinf(x); }
inline float __expf(float x) { return expf(x); }
inline float __logf(float x) { return logf(x); }
inline float __powf(float x, float y) { return powf(x, y); }
inline float __fdividef(float x, float y) { return x / y; }
inline float rsqrtf(float x) { return 1.f / sqrtf(x); }
inline double rsqrt(double x) { return 1. / sqrt(x); }
#endif
#ifdef HOST_EMULATION
#include <cstring>
inline int __float_as_int(float f) { int i; memcpy(&i, &f, 4); return i; }
inline float __int_as_float(int i) { float f; memcpy(&f, &i, 4); return f; }
#endif

#ifdef HOST_EMULATION
#include <cstdlib>
template<class T> cudaError_t cudaMemcpyToSymbol(const T &sym, const void *src, size_t size, size_t offset = 0, cudaMemcpyKind = cudaMemcpyHostToDevice) { memcpy((char *) &sym + offset, src, size); return 0; }
inline cudaError_t cudaMemcpy(void *dst, const void *src, size_t size, cudaMemcpyKind) { memcpy(dst, src, size); return 0; }
inline cudaError_t cudaMalloc(void **p, size_t s) { *p = calloc(1, s); return 0; }
inline cudaError_t cudaFree(void *p) { free(p); return 0; }
template<class T> cudaError_t cudaMalloc(T **p, size_t s) { return cudaMalloc((void**)p, s); }
#endif
