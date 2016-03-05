#ifndef FDTD_H
#define FDTD_H

#include <cuda_runtime_api.h>
#include <pthread.h>

typedef struct
{
    int   nx, ny, nz;
    int   iterationsCount;
    char  *inputPath;
    char  *outputPath;
    int   elementsPerWave;
    float waveFrequency;
    float pulseWidth;
    float pulseModulationFrequency;
    int   sourcesCount;
    int   *sources;
    float defaultSigma;
    float defaultEpsS;
    float defaultEpsI;
    float defaultTauD;
    int   jzCount;
    float *jz;

    float pi; // Delicious pie
    float c; // Light speed (v in m/s)
    float timeskip; // Time step skip
    float lambda; // Wave length (meters)
    float dt; // Length of the time step
    float dx, dy, dz; // Distance between 2 cells
    float mu0; // Permeability of free space (in henry/meter)
    float eps0; // Permittivity of free space (in farad/meter)
} FdtdParams;


typedef struct
{
    float *ex0, *ey0, *ez0;
    float *ex1, *ey1, *ez1;
    float *ex2, *ey2, *ez2;

    float *hx, *hy, *hz;

    float *dx0, *dy0, *dz0;
    float *dx1, *dy1, *dz1;
    float *dx2, *dy2, *dz2;

    float *epsI, *epsS;
    float *tauD, *sigma;

    float *rpx0,   *rpy0,   *rpz0;
    float *rpxEnd, *rpyEnd, *rpzEnd;
} FdtdField;


typedef struct
{
    float *xSource;
    float *ySource;
    float *zSource;

    float *xBuffer;
    float *yBuffer;
    float *zBuffer;

    FdtdParams *params;
    cudaStream_t stream;
    pthread_mutex_t *copyMutex;
} CopyParams;


typedef struct
{
    FdtdParams *params;

    CopyParams *hParams;
    CopyParams *dParams;
    CopyParams *eParams;

    int currentIteration;
} ResultsParams;


FdtdParams *initParamsWithPath(const char *);
void deallocParams(FdtdParams *);
void printParams(FdtdParams *params);

FdtdField *initFieldWithParams(FdtdParams *);
void deallocField(FdtdField *);

FdtdField *initDeviceFieldWithParams(FdtdParams *);
void deallocDeviceField(FdtdField *);

void setupMurBoundary(FdtdParams *params, FdtdField *field);
void loadMaterials(FdtdParams *params, FdtdField *field, const char *specsFilePath, const char *materialsPath);
void setupSources(FdtdParams *params);
void copyDataToDevice(FdtdParams *params, FdtdField *field, FdtdField *deviceField);
void *copyResultsWithParams(void *params);
void *writeResultsWithParams(void *params);
void writeResults(FdtdParams *params,
                  float *hxSource, float *hySource, float *hzSource,
                  float *dxSource, float *dySource, float *dzSource,
                  float *exSource, float *eySource, float *ezSource,
                  int currentIteration);
void writeResults3d(FdtdParams *params,
                    float *hxSource, float *hySource, float *hzSource,
                    float *dxSource, float *dySource, float *dzSource,
                    float *exSource, float *eySource, float *ezSource,
                    int currentIteration);

#endif
