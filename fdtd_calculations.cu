#include "fdtd_calculations.h"


#define OFFSET(p, x, y, z) (p[(z)*ny*nx + (y)*nx + (x)])


__global__ void updateHField(float *hx,       float *hy,       float *hz,
                             float *exSource, float *eySource, float *ezSource,
                             int nx, int ny, int nz,
                             float dt, float dx, float dy, float dz,
                             float mu0)
{
    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int iz = threadIdx.z + blockIdx.z * blockDim.z;

    // Update hx
    if(ix > 0 && ix < nx-1 &&
       iy > 0 && iy < ny-1 &&
       iz > 0 && iz < nz-1) {
        OFFSET(hx, ix, iy, iz) = OFFSET(hx, ix, iy, iz) -
                                 dt/(mu0 * dy) *
                                 (OFFSET(ezSource, ix, iy+1, iz) - OFFSET(ezSource, ix, iy, iz)) +
                                 dt/(mu0 * dz) *
                                 (OFFSET(eySource, ix, iy, iz+1) - OFFSET(eySource, ix, iy, iz));
    }
    
    // Update hy
    if(ix > 0 && ix < nx-1 &&
       iy > 1 && iy < ny-1 &&
       iz > 0 && iz < nz-1) {
        OFFSET(hy, ix, iy, iz) = OFFSET(hy, ix, iy, iz) -
                                 dt/(mu0 * dz) *
                                 (OFFSET(exSource, ix, iy, iz+1) - OFFSET(exSource, ix, iy, iz)) +
                                 dt/(mu0 * dx) *
                                 (OFFSET(ezSource, ix+1, iy, iz) - OFFSET(ezSource, ix, iy, iz));
    }
    
    // Update hz
    if(ix > 0 && ix < nx-1 &&
       iy > 0 && iy < ny-1 &&
       iz > 1 && iz < nz-1) {
        OFFSET(hz, ix, iy, iz) = OFFSET(hz, ix, iy, iz) -
                                 dt/(mu0 * dx) *
                                 (OFFSET(eySource, ix+1, iy, iz) - OFFSET(eySource, ix, iy, iz)) +
                                 dt/(mu0 * dy) *
                                 (OFFSET(exSource, ix, iy+1, iz) - OFFSET(exSource, ix, iy, iz));
    }
}


__global__ void updateDField(float *dxTarget, float *dyTarget, float *dzTarget,
                             float *dxSource, float *dySource, float *dzSource,
                             float *hx,       float *hy,       float *hz,
                             int nx, int ny, int nz,
                             float dt, float dx, float dy, float dz)
{
    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int iz = threadIdx.z + blockIdx.z * blockDim.z;

    // Update dx
    if(ix > 0 && ix < nx-1 &&
       iy > 1 && iy < ny-1 &&
       iz > 1 && iz < nz-1) {
        OFFSET(dxSource, ix, iy, iz) = OFFSET(dxSource, ix, iy, iz) +
                                       dt/dy * (OFFSET(hz, ix, iy, iz) - OFFSET(hz, ix, iy-1, iz)) -
                                       dt/dz * (OFFSET(hy, ix, iy, iz) - OFFSET(hy, ix, iy, iz-1));
    }
    
    // Update dy
    if(ix > 1 && ix < nx-1 &&
       iy > 0 && iy < ny-1 &&
       iz > 1 && iz < nz-1) {
        OFFSET(dyTarget, ix, iy, iz) = OFFSET(dySource, ix, iy, iz) +
                                       dt/dz * (OFFSET(hx, ix, iy, iz) - OFFSET(hx, ix, iy, iz-1)) -
                                       dt/dx * (OFFSET(hz, ix, iy, iz) - OFFSET(hz, ix-1, iy, iz));
    }
    
    // Update dz
    if(ix > 1 && ix < nx-1 &&
       iy > 1 && iy < ny-1 &&
       iz > 0 && iz < nz-1) {
            OFFSET(dzTarget, ix, iy, iz) = OFFSET(dzSource, ix, iy, iz) +
                                           dt/dx * (OFFSET(hy, ix, iy, iz) - OFFSET(hy, ix-1, iy, iz)) -
                                           dt/dy * (OFFSET(hx, ix, iy, iz) - OFFSET(hx, ix, iy-1, iz));
    }
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
    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int iz = threadIdx.z + blockIdx.z * blockDim.z;

    // Update ex
    if(ix > 0 && ix < nx-1 &&
       iy > 1 && iy < ny-1 &&
       iz > 1 && iz < nz-1) {
        OFFSET(exTarget, ix, iy, iz) = (
                                        1.0/(2.0 * eps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz) +
                                        2.0 * dt *
                                        (
                                         eps0 * OFFSET(epsS, ix, iy, iz) +
                                         OFFSET(sigma, ix, iy, iz) * OFFSET(tauD, ix, iy, iz)
                                        ) +
                                        OFFSET(sigma, ix, iy, iz) * dt * dt)
                                       ) *
                                       (
                                        (
                                         4.0 * eps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz) +
                                         2.0 * dt *
                                         (
                                          eps0 * OFFSET(epsS, ix, iy, iz) +
                                          OFFSET(sigma, ix, iy, iz) * OFFSET(tauD, ix, iy, iz)
                                         ) -
                                         OFFSET(sigma, ix, iy, iz) * dt * dt
                                        ) *
                                        OFFSET(exSource0, ix, iy, iz) -
                                        (2.0 * eps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz)) *
                                        OFFSET(exSource1, ix, iy, iz) +
                                        (2.0 * (dt + OFFSET(tauD, ix, iy, iz))) * OFFSET(dxSource0, ix, iy, iz) -
                                        (2.0 * dt + 4.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dxSource1, ix, iy, iz) +
                                        (2.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dxSource2, ix, iy, iz)
                                       );
    }
    
    // Update ey
    if(ix >= 2 && ix <= nx-1 &&
       iy >= 1 && iy <= ny-1 &&
       iz >= 2 && iz <= nz-1) {
        OFFSET(eyTarget, ix, iy, iz) = (
                                        1.0/(2.0 * eps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz) +
                                        2.0 * dt *
                                        (
                                         eps0 * OFFSET(epsS, ix, iy, iz) +
                                         OFFSET(sigma, ix, iy, iz) * OFFSET(tauD, ix, iy, iz)
                                        ) +
                                        OFFSET(sigma, ix, iy, iz) * dt * dt)
                                       ) *
                                       (
                                        (
                                         4.0 * eps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz) +
                                         2.0 * dt *
                                         (
                                          eps0 * OFFSET(epsS, ix, iy, iz) +
                                          OFFSET(sigma, ix, iy, iz) * OFFSET(tauD, ix, iy, iz)
                                         ) -
                                         OFFSET(sigma, ix, iy, iz) * dt * dt
                                        ) *
                                        OFFSET(eySource0, ix, iy, iz) -
                                        (2.0 * eps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz)) *
                                        OFFSET(eySource1, ix, iy, iz) +
                                        (2.0 * (dt + OFFSET(tauD, ix, iy, iz))) * OFFSET(dySource0, ix, iy, iz) -
                                        (2.0 * dt + 4.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dySource1, ix, iy, iz) +
                                        (2.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dySource2, ix, iy, iz)
                                       );
    }
    
    // Update ez
    if(ix >= 2 && ix <= nx-1 &&
       iy >= 1 && iy <= ny-1 &&
       iz >= 2 && iz <= nz-1) {
        OFFSET(ezTarget, ix, iy, iz) = (
                                        1.0/(2.0 * eps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz) +
                                        2.0 * dt *
                                        (
                                         eps0 * OFFSET(epsS, ix, iy, iz) +
                                         OFFSET(sigma, ix, iy, iz) * OFFSET(tauD, ix, iy, iz)
                                        ) +
                                        OFFSET(sigma, ix, iy, iz) * dt * dt)
                                       ) *
                                       (
                                        (
                                         4.0 * eps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz) +
                                         2.0 * dt *
                                         (
                                          eps0 * OFFSET(epsS, ix, iy, iz) +
                                          OFFSET(sigma, ix, iy, iz) * OFFSET(tauD, ix, iy, iz)
                                         ) -
                                         OFFSET(sigma, ix, iy, iz) * dt * dt
                                        ) *
                                        OFFSET(ezSource0, ix, iy, iz) -
                                        (2.0 * eps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz)) *
                                        OFFSET(ezSource1, ix, iy, iz) +
                                        (2.0 * (dt + OFFSET(tauD, ix, iy, iz))) * OFFSET(dzSource0, ix, iy, iz) -
                                        (2.0 * dt + 4.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dzSource1, ix, iy, iz) +
                                        (2.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dzSource2, ix, iy, iz)
                                       );
    }
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
