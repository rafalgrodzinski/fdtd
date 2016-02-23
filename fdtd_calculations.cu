#include "fdtd_calculations.h"

#include <stdio.h>
#include "utils.h"


void copySymbolsToDevice(FdtdParams *params)
{

    CHECK(cudaMemcpyToSymbol(dNx, &params->nx, sizeof(int)));
    CHECK(cudaMemcpyToSymbol(dNy, &params->ny, sizeof(int)));
    CHECK(cudaMemcpyToSymbol(dNz, &params->nz, sizeof(int)));

    CHECK(cudaMemcpyToSymbol(dDt, &params->dt, sizeof(float)));
    CHECK(cudaMemcpyToSymbol(dDx, &params->dx, sizeof(float)));
    CHECK(cudaMemcpyToSymbol(dDy, &params->dy, sizeof(float)));
    CHECK(cudaMemcpyToSymbol(dDz, &params->dz, sizeof(float)));

    CHECK(cudaMemcpyToSymbol(dMu0, &params->mu0, sizeof(float)));
    CHECK(cudaMemcpyToSymbol(dEps0, &params->eps0, sizeof(float)));

    CHECK(cudaMemcpyToSymbol(dSourcesCount, &params->sourcesCount, sizeof(int)));

    int sourcesCount = (D_MAX_SOURCES < params->sourcesCount) ? D_MAX_SOURCES : params->sourcesCount;
    CHECK(cudaMemcpyToSymbol(dSources, params->sources, sizeof(int) * 3 * sourcesCount));

    int jzCount = (D_MAX_JZ < params->jzCount) ? D_MAX_JZ : params->jzCount;
    CHECK(cudaMemcpyToSymbol(dJz, params->jz, sizeof(float) * jzCount));
}


__global__ void updateHField(float *hx,       float *hy,       float *hz,
                             float *exSource, float *eySource, float *ezSource)
{
    int nx = dNx;
    int ny = dNy;
    int nz = dNz;

    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int iz = threadIdx.z + blockIdx.z * blockDim.z;

    // Update hx
    if(ix >= 1 && ix < nx-1 &&
       iy >= 0 && iy < ny-1 &&
       iz >= 0 && iz < nz-1) {

        OFFSET(hx, ix, iy, iz) = OFFSET(hx, ix, iy, iz) -
                                 dDt/(dMu0 * dDy) *
                                 (OFFSET(ezSource, ix, iy+1, iz) - OFFSET(ezSource, ix, iy, iz)) +
                                 dDt/(dMu0 * dDy) *
                                 (OFFSET(eySource, ix, iy, iz+1) - OFFSET(eySource, ix, iy, iz));
    }

    
    // Update hy
    if(ix >= 0 && ix < nx-1 &&
       iy >= 1 && iy < ny-1 &&
       iz >= 0 && iz < nz-1) {
        OFFSET(hy, ix, iy, iz) = OFFSET(hy, ix, iy, iz) -
                                 dDt/(dMu0 * dDz) *
                                 (OFFSET(exSource, ix, iy, iz+1) - OFFSET(exSource, ix, iy, iz)) +
                                 dDt/(dMu0 * dDx) *
                                 (OFFSET(ezSource, ix+1, iy, iz) - OFFSET(ezSource, ix, iy, iz));
    }
    
    // Update hz
    if(ix >= 0 && ix < nx-1 &&
       iy >= 0 && iy < ny-1 &&
       iz >= 1 && iz < nz-1) {
        OFFSET(hz, ix, iy, iz) = OFFSET(hz, ix, iy, iz) -
                                 dDt/(dMu0 * dDx) *
                                 (OFFSET(eySource, ix+1, iy, iz) - OFFSET(eySource, ix, iy, iz)) +
                                 dDt/(dMu0 * dDy) *
                                 (OFFSET(exSource, ix, iy+1, iz) - OFFSET(exSource, ix, iy, iz));
    }
}


__global__ void updateDField(float *dxTarget, float *dyTarget, float *dzTarget,
                             float *dxSource, float *dySource, float *dzSource,
                             float *hx,       float *hy,       float *hz)
{
    int nx = dNx;
    int ny = dNy;
    int nz = dNz;

    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int iz = threadIdx.z + blockIdx.z * blockDim.z;

    float x = dDt/dDx;
    float y = dDt/dDy;
    float z = dDt/dDz;

    // Update dDx
    if(ix >= 0 && ix < nx-1 &&
       iy >= 1 && iy < ny-1 &&
       iz >= 1 && iz < nz-1) {
        OFFSET(dxTarget, ix, iy, iz) = OFFSET(dxSource, ix, iy, iz) +
                                       y * (OFFSET(hz, ix, iy, iz) - OFFSET(hz, ix, iy-1, iz)) -
                                       z * (OFFSET(hy, ix, iy, iz) - OFFSET(hy, ix, iy, iz-1));
    }
    
    // Update dDy
    if(ix >= 1 && ix < nx-1 &&
       iy >= 0 && iy < ny-1 &&
       iz >= 1 && iz < nz-1) {
        OFFSET(dyTarget, ix, iy, iz) = OFFSET(dySource, ix, iy, iz) +
                                       z * (OFFSET(hx, ix, iy, iz) - OFFSET(hx, ix, iy, iz-1)) -
                                       x * (OFFSET(hz, ix, iy, iz) - OFFSET(hz, ix-1, iy, iz));
    }
    
    // Update dDz
    if(ix >= 1 && ix < nx-1 &&
       iy >= 1 && iy < ny-1 &&
       iz >= 0 && iz < nz-1) {
            OFFSET(dzTarget, ix, iy, iz) = OFFSET(dzSource, ix, iy, iz) +
                                           x * (OFFSET(hy, ix, iy, iz) - OFFSET(hy, ix-1, iy, iz)) -
                                           y * (OFFSET(hx, ix, iy, iz) - OFFSET(hx, ix, iy-1, iz));
    }
}


__global__ void updateEField(float *exTarget,  float *eyTarget,  float *ezTarget,
                             float *exSource0, float *eySource0, float *ezSource0,
                             float *exSource1, float *eySource1, float *ezSource1,
                             float *dxSource0, float *dySource0, float *dzSource0,
                             float *dxSource1, float *dySource1, float *dzSource1,
                             float *dxSource2, float *dySource2, float *dzSource2,
                             float *sigma, float *epsI, float *epsS, float *tauD)
{
    int nx = dNx;
    int ny = dNy;
    int nz = dNz;

    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int iz = threadIdx.z + blockIdx.z * blockDim.z;

    float a = 2.0 * dEps0 * OFFSET(epsI, ix, iy, iz) * OFFSET(tauD, ix, iy, iz);
    float b = a + 2.0 * dDt * (dEps0 * OFFSET(epsS, ix, iy, iz) + OFFSET(sigma, ix, iy, iz) * OFFSET(tauD, ix, iy, iz));
    float c = 1.0/(b + OFFSET(sigma, ix, iy, iz) * dDt * dDt);
    float d = (2.0 * b - OFFSET(sigma, ix, iy, iz) * dDt * dDt);

    // Update ex
    if(ix >= 0 && ix < nx-1 &&
       iy >= 1 && iy < ny-1 &&
       iz >= 1 && iz < nz-1) {
        OFFSET(exTarget, ix, iy, iz) = c *
                                       (d * OFFSET(exSource0, ix, iy, iz) -
                                        a * OFFSET(exSource1, ix, iy, iz) +
                                        (2.0 * (dDt + OFFSET(tauD, ix, iy, iz))) * OFFSET(dxSource0, ix, iy, iz) -
                                        (2.0 * dDt + 4.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dxSource1, ix, iy, iz) +
                                        (2.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dxSource2, ix, iy, iz)
                                       );
    }
    
    // Update ey
    if(ix >= 1 && ix <= nx-1 &&
       iy >= 0 && iy <= ny-1 &&
       iz >= 1 && iz <= nz-1) {
        OFFSET(eyTarget, ix, iy, iz) = c *
                                       (d * OFFSET(eySource0, ix, iy, iz) -
                                        a * OFFSET(eySource1, ix, iy, iz) +
                                        (2.0 * (dDt + OFFSET(tauD, ix, iy, iz))) * OFFSET(dySource0, ix, iy, iz) -
                                        (2.0 * dDt + 4.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dySource1, ix, iy, iz) +
                                        (2.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dySource2, ix, iy, iz)
                                       );
    }
    
    // Update ez
    if(ix >= 1 && ix <= nx-1 &&
       iy >= 0 && iy <= ny-1 &&
       iz >= 1 && iz <= nz-1) {
        OFFSET(ezTarget, ix, iy, iz) = c *
                                       (d * OFFSET(ezSource0, ix, iy, iz) -
                                        a * OFFSET(ezSource1, ix, iy, iz) +
                                        (2.0 * (dDt + OFFSET(tauD, ix, iy, iz))) * OFFSET(dzSource0, ix, iy, iz) -
                                        (2.0 * dDt + 4.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dzSource1, ix, iy, iz) +
                                        (2.0 * OFFSET(tauD, ix, iy, iz)) * OFFSET(dzSource2, ix, iy, iz)
                                       );
    }
}


__global__ void updateSources(float *dzTarget, float *dzSource,
                              float *hx,       float *hy,
                              int currIteration)
{
    int nx = dNx;
    int ny = dNy;

    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int iz = threadIdx.z + blockIdx.z * blockDim.z;

    // Update source
    if(ix == 0 && iy == 0 && iz == 0) {
        for(int i=0; i < dSourcesCount; i++) {
            int x = dSources[i * 3 + 0];
            int y = dSources[i * 3 + 1];
            int z = dSources[i * 3 + 2];

            float val = OFFSET(dzSource, x, y, z) +
                                        dDt/dDx * (OFFSET(hy, x, y, z) - OFFSET(hy, x-1, y, z)) -
                                        dDt/dDy * (OFFSET(hx, x, y, z) - OFFSET(hx, x, y-1, z)) -
                                        dJz[currIteration];

            OFFSET(dzTarget, x, y, z) =  val;        
        }
    }
}


__global__ void updateMurBoundary(float *exTarget, float *eyTarget, float *ezTarget,
                                  float *exSource, float *eySource, float *ezSource,
                                  float *rpx0,     float *rpy0,     float *rpz0,
                                  float *rpxEnd,   float *rpyEnd,   float *rpzEnd)
{
    int nx = dNx;
    int ny = dNy;
    int nz = dNz;

    int rpnx, rpny;

    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int iz = threadIdx.z + blockIdx.z * blockDim.z;
    
    // Update ex

    // for rpy
    rpnx = dNx;
    rpny = 2;

    if(ix >= 0 && ix < nx-1 && 
       iy == 0 && 
       iz >= 1 && iz < nz-1) {
        OFFSET(exTarget, ix, iy, iz) = 1/(dDt + dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpy0, ix, 0, iz))) *  
                                       (                                                        
                                        (dDt - dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpy0, ix, 1, iz))) * 
                                        OFFSET(exTarget, ix, iy+1, iz) +                               
                                        (dDt + dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpy0, ix, 1, iz))) * 
                                        OFFSET(exSource, ix, iy+1, iz) -                               
                                        (dDt - dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpy0, ix, 0, iz))) *    
                                        OFFSET(exSource, ix, iy, iz)                                   
                                       );
    }

    if(ix >= 0    && ix < nx-1 && 
       iy == ny-1 && 
       iz >= 1    && iz < nz-1) {
        OFFSET(exTarget, ix, iy, iz) = 1/(dDt + dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpyEnd, ix, 1, iz))) *  
                                       (                                                          
                                        (dDt - dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpyEnd, ix, 0, iz))) * 
                                        OFFSET(exTarget, ix, iy-1, iz) +                                 
                                        (dDt + dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpyEnd, ix, 0, iz))) * 
                                        OFFSET(exSource, ix, iy-1, iz) -                                 
                                        (dDt - dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpyEnd, ix, 1, iz))) *   
                                        OFFSET(exSource, ix, iy, iz)                                     
                                       );
    }

    // for rpz
    rpnx = dNx;
    rpny = dNy;

    if(ix >= 0 && ix < nx-1 && 
       iy >= 1 && iy < ny-1 && 
       iz == 0) {
        OFFSET(exTarget, ix, iy, iz) = 1/(dDt + dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpz0, ix, iy, 0))) *  
                                       (                                                        
                                        (dDt - dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpz0, ix, iy, 1))) * 
                                        OFFSET(exTarget, ix, iy, iz+1) +                               
                                        (dDt + dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpz0, ix, iy, 1))) * 
                                        OFFSET(exSource, ix, iy, iz+1) -                               
                                        (dDt - dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpz0, ix, iy, 0))) *   
                                        OFFSET(exSource, ix, iy, iz)                                   
                                       );
    }

    if(ix >= 0 && ix < nx-1 && 
       iy >= 1 && iy < ny-1 && 
       iz == nz-1) {
        OFFSET(exTarget, ix, iy, iz) = 1/(dDt + dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpzEnd, ix, iy, 1))) *  
                                       (                                                          
                                        (dDt - dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpzEnd, ix, iy, 0))) * 
                                        OFFSET(exTarget, ix, iy, iz-1) +                                 
                                        (dDt + dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpzEnd, ix, iy, 0))) * 
                                        OFFSET(exSource, ix, iy, iz-1) -                                 
                                        (dDt - dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpzEnd, ix, iy, 1))) *   
                                        OFFSET(exSource, ix, iy, iz)                                     
                                       );
    }

    // Update ey

    // for rpx
    rpnx = 2;
    rpny = dNy;

    if(ix == 0 && 
       iy >= 0  && iy < ny-1 && 
       iz >= 1  && iz < nz-1) {
        OFFSET(eyTarget, ix, iy, iz) = 1/(dDt + dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpx0, 0, iy, iz))) *  
                                       (                                                        
                                        (dDt - dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpx0, 1, iy, iz))) * 
                                        OFFSET(eyTarget, ix+1, iy, iz) +                               
                                        (dDt + dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpx0, 1, iy, iz))) * 
                                        OFFSET(eySource, ix+1, iy, iz) -                               
                                        (dDt - dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpx0, 0, iy, iz))) *   
                                        OFFSET(eySource, ix, iy, iz)                                   
                                       );
    }

    if(ix == nx-1 && 
       iy >= 0     && iy < ny-1 && 
       iz >= 1     && iz < nz-1) {
        OFFSET(eyTarget, ix, iy, iz) = 1/(dDt + dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpxEnd, 1, iy, iz)))  * 
                                       (                                                          
                                        (dDt - dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpxEnd, 0, iy, iz))) * 
                                        OFFSET(eySource, ix-1, iy, iz) +                                 
                                        (dDt + dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpxEnd, 0, iy, iz))) * 
                                        OFFSET(eySource, ix-1, iy, iz) -                                 
                                        (dDt - dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpxEnd, 1, iy, iz))) *   
                                        OFFSET(eySource, ix, iy, iz)                                     
                                       );
    }

    // for rpz
    rpnx = dNx;
    rpny = dNy;

    if(ix >= 1 && ix < nx-1 && 
       iy >= 0 && iy < ny-1 && 
       iz == 0) {
        OFFSET(eyTarget, ix, iy, iz) = 1/(dDt + dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpz0, ix, iy, 0))) *  
                                       (                                                        
                                        (dDt - dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpz0, ix, iy, 1))) * 
                                        OFFSET(eyTarget, ix, iy,iz+1) +                                
                                        (dDt + dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpz0, ix, iy, 1))) * 
                                        OFFSET(eySource, ix, iy, iz+1) -                               
                                        (dDt - dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpz0, ix, iy, 0))) *   
                                        OFFSET(eySource, ix, iy, iz)                                   
                                       );
    }

    if(ix >= 1 && ix < nx-1 && 
       iy >= 0 && iy < ny-1 && 
       iz == nz-1) {
        OFFSET(eyTarget, ix, iy, iz) = 1/(dDt + dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpzEnd, ix, iy, 1))) *  
                                       (                                                          
                                        (dDt - dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpzEnd, ix, iy, 0))) * 
                                        OFFSET(eyTarget, ix, iy, iz-1) +                                 
                                        (dDt + dDz * sqrt(dMu0 * dEps0 * OFFSETRP(rpzEnd, ix, iy, 0))) * 
                                        OFFSET(eySource, ix, iy, iz-1) -                                 
                                        (dDt - dDz *sqrt(dMu0 * dEps0 * OFFSETRP(rpzEnd, ix, iy, 1))) *    
                                        OFFSET(eySource, ix, iy, iz)                                     
                                       );
    }

    // Update ez

    // for rpz
    rpnx = 2;
    rpny = dNy;

    if(ix == 0 && 
       iy >= 1  && iy < ny-1 && 
       iz >= 0  && iz < nz-1) {
        OFFSET(ezTarget, ix, iy, iz) = 1/(dDt + dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpx0, 0, iy, iz))) *  
                                       (                                                        
                                        (dDt - dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpx0, 1, iy, iz))) * 
                                        OFFSET(ezTarget, ix+1, iy, iz) +                               
                                        (dDt + dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpx0, 1, iy, iz))) * 
                                        OFFSET(ezSource, ix+1, iy, iz) -                                
                                        (dDt - dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpx0, 0, iy, iz)))  *  
                                        OFFSET(ezSource, ix, iy, iz)                                   
                                       );
    }
      
    if(ix == nx-1 && 
       iy >= 1     && iy < ny-1 && 
       iz >= 0     && iz < nz-1) {
        OFFSET(ezTarget, ix, iy, iz) = 1/(dDt + dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpxEnd, 1, iy, iz))) *  
                                       (                                                          
                                        (dDt - dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpxEnd, 0, iy, iz))) * 
                                        OFFSET(ezTarget, ix-1, iy, iz) +                                 
                                        (dDt + dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpxEnd, 0, iy, iz))) * 
                                        OFFSET(ezSource, ix-1, iy, iz) -                                 
                                        (dDt - dDx * sqrt(dMu0 * dEps0 * OFFSETRP(rpxEnd, 1, iy, iz))) *   
                                        OFFSET(ezSource, ix, iy, iz)                                     
                                       );
    }

    // for rpy
    rpnx = dNx;
    rpny = 2;
    
    if(ix >= 1  && ix < nx-1 && 
       iy == 0 && 
       iz >= 0  && iz < nz-1) { 
        OFFSET(ezTarget, ix, iy, iz) = 1/(dDt + dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpy0, ix, 0, iz))) *  
                                       (                                                        
                                        (dDt - dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpy0, ix, 1, iz))) * 
                                        OFFSET(ezTarget, ix, iy+1, iz) +                                
                                        (dDt + dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpy0, ix, 1, iz))) * 
                                        OFFSET(ezSource, ix, iy+1, iz) -                               
                                        (dDt - dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpy0, ix, 0, iz))) *   
                                        OFFSET(ezSource, ix, iy, iz)                                   
                                       );
    }
      
    if(ix >= 1     && ix < nx-1 && 
       iy == ny-1 && 
       iz >= 0     && iz < nz-1) { 
            OFFSET(ezTarget, ix, iy, iz) = 1/(dDt + dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpyEnd, ix, 1, iz))) *  
                                           (                                                          
                                            (dDt - dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpyEnd, ix, 0, iz))) * 
                                            OFFSET(ezTarget, ix, iy-1, iz) +                                 
                                            (dDt + dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpyEnd, ix, 0, iz))) * 
                                            OFFSET(ezSource, ix, iy-1, iz) -                                 
                                            (dDt - dDy * sqrt(dMu0 * dEps0 * OFFSETRP(rpyEnd, ix, 1, iz))) *   
                                            OFFSET(ezSource, ix, iy, iz)                                     
                                           );
    }
}
