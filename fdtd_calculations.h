#ifndef FDTD_CALCULATIONS_H
#define FDTD_CALCULATIONS_H

#include <cuda_runtime_api.h>
#include "fdtd.h"

#define D_MAX_SOURCES 16
#define D_MAX_JZ 512

__constant__ int dNx, dNy, dNz;
__constant__ float dDt, dDx, dDy, dDz;
__constant__ float dMu0, dEps0;
__constant__ int dSourcesCount;
__constant__ float dSources[D_MAX_SOURCES * 3];
__constant__ float dJz[D_MAX_JZ];


void copySymbolsToDevice(FdtdParams *params);

__global__ void updateHField(float *hx,       float *hy,       float *hz,
                             float *exSource, float *eySource, float *ezSource);

__global__ void updateDField(float *dxTarget, float *dyTarget, float *dzTarget,
                             float *dxSource, float *dySource, float *dzSource,
                             float *hx,       float *hy,       float *hz);

__global__ void updateEField(float *exTarget,  float *eyTarget,  float *ezTarget,
                             float *exSource0, float *eySource0, float *ezSource0,
                             float *exSource1, float *eySource1, float *ezSource1,
                             float *dxSource0, float *dySource0, float *dzSource0,
                             float *dxSource1, float *dySource1, float *dzSource1,
                             float *dxSource2, float *dySource2, float *dzSource2,
                             float *sigma,float *epsI, float *epsS, float *tauD);

__global__ void updateSources(float *dzTarget, float *dzSource,
                              float *hx,       float *hy,
                              int currIteration);

__global__ void updateMurBoundary(float *exTarget, float *eyTarget, float *ezTarget,
                                  float *exSource, float *eySource, float *ezSource,
                                  float *rpx0,     float *rpy0,     float *rpz0,
                                  float *rpxEnd,   float *rpyEnd,   float *rpZend);

#endif
