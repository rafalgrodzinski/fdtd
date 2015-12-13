#include "fdtd_calculations.h"

__global__ void updateHField(float *hx,       float *hy,       float *hz,
                             float *exSource, float *eySource, float *ezSource,
                             int nx, int ny, int nz,
                             float dt, float dx, float dy, float dz,
                             float mu0)
{
}


__global__ void updateDField(float *dxTarget, float *dyTarget, float *dzTarget,
                             float *dxSource, float *dySource, float *dzSource,
                             float *hx,       float *hy,       float *hz,
                             int nx, int ny, int nz,
                             float dt, float dx, float dy, float dz)
{
}


__global__ void updateEField(float *exTarget,  float *eyTarget,  float *ezTarget,
                             float *exSource0, float *eySource0, float *ezSource0,
                             float *exSource1, float *eySource1, float *ezSource1,
                             float *dxSource0, float *dySource0, float *dzSource0,
                             float *dxSource1, float *dySource1, float *dzSource1,
                             float *dxSource2, float *dySource2, float *dzSource2,
                             float *sigma,float *epsI, float *epsS, float *tauD,
                             int nx, int ny, int nz,
                             float dt, float eps0)
{
}


__global__ void updateSource(float *dzTarget, float *dzSource,
                             float *hx,       float *hy,
                             int *src, float *jz,
                             float dt, float dx, float dy, float dz,
                             int nsrc, int runsCount)
{
}


__global__ void updateMurBoundary(float *exTarget, float *eyTarget, float *ezTarget,
                                  float *exSource, float *eySource, float *ezSource,
                                  float *rpx0,     float *rpy0,     float *rpz0,
                                  float *rpxEnd,   float *rpyEnd,   float *rpZend,
                                  int nx, int ny, int nz,
                                  float dt, float dx, float dy, float dz,
                                  float mu0, float eps0)
{
}
