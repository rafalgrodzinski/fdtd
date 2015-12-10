#ifndef FDTD_H
#define FDTD_H

typedef struct
{
    int nx, ny, nz;
    int iterationsCount;
    char *inputPath;
    char *outputPath;
    int  elementsPerWave;
    float waveFrequency;
    float pulseWidth;
    float pulseModulationFrequency;
    int   sourcesCount;
    int   *sources;
    float sigma;
    float eps_s;
    float eps_i;
    float tau_d;
} FdtdParams;


typedef struct
{
    float *ex1, *ey1, *ez1;
    float *ex2, *ey2, *ez2;
    float *ex3, *ey3, *ez3;

    float *hx, *hy, *hz;

    float *dx1, *dy1, *dz1;
    float *dx2, *dy2, *dz2;
    float *dx3, *dy3, *dz3;

    float *eps_i, *eps_s;
    float *tau_d, *sigma;

    float *rp_x_0,   *rp_y_0,   *rp_z_0;
    float *rp_x_end, *rp_y_end, *rp_z_end;
} FdtdField;


void printUsage();
FdtdParams *initParamsWithPath(const char *);
void deallocParams(FdtdParams *);

FdtdField  *initHostFieldWithParams(FdtdParams *);
void deallocHostField(FdtdField *);

FdtdField  *initDeviceFieldWithParams(FdtdParams *);
void deallocDeviceField(FdtdField *);

void printParams(FdtdParams *params);

#endif
