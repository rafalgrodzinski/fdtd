#include <stdio.h>


#define BLOCK_X 4
#define BLOCK_Y 4
#define BLOKC_Z 4


typedef struct FdtdParams
{
    int nx, ny, nz;
    int iterationsCount;
    char *inputPath;
    char *outputPath;
    int  elementsPerWave;
    float waveFrequency;
    float pulseWidth;
    float pulseModulationFrequency;
    int   sourceCount;
    int   *sources;
    float sigma;
    float eps_i;
    float tau_d;
};


typedef struct FdtdField
{
    float *ex1, *ey1, *ez1;
    float *ex2, *ey2, *ez2;
    float *ex3, *ey3, *ez3;

    float *hx, *hy, *hz;

    float *dx1, *dy1, *dz1;
    float *dx2, *dy2, *dz2;
    float *dx3, *dy3, *dz3;

    float *eps_i, *eps_s;
    float *tau_d, *sigma

    float *rp_x_0,   *rp_y_0,   *rp_z_0;
    float *rp_x_end, *rp_y_end, *rp_z_end;
};


void printUsage();
FdtdParams *initParamsWithPath(char *);

FdtdField  *initHostFieldWithParams(FdtdParams *);
void deallocHostField(FdtdField *);

FdtdField  *initDeviceFieldWithParams(FdtdParams *);
void deallocDeviceField(FdtdField *);


int main(int argc, char **argv)
{
    // Read params and initialize field
    FdtdParams *params;
    FdtdField  *hostField, *deviceField; // Used for CUDA

    params = initParamsWithPath("data/input_params");

    hostField = initHostFieldWithParams(params);
    deviceField = initDeviceFieldWithParams(params);

    // Setup CUDA parameters
    dim3 gridSize((params->nx + BLOCK_X - 1)/BLOCK_X,
                  (params->ny + BLOCK_Y - 1)/BLOCK_Y,
                  (params->nz + BLOCK_Z - 1)/BLOCK_Z);
    dim3 blockSize(BLOCK_X, BLOCK_Y, BLOCK_Z);

    // Main loop
    for(int i=0; i<params->iterationsCount; i++) {
        // Run 0
        // Run 1
        // Run 2
    }
}


void printUsage()
{
}


FdtdParams *initParamsWithPath(char *filePath)
{
    FdtdParams *params = malloc(sizeof(FdtdParams));

    FILE paramsFile = fopen(filePath, "r");
    
    char temp[100];

    // nx, ny, nz
    fscanf(paramsFile, "%s %d %d %d", temp, &params->nx, &params->ny, &params->nz);
    // runsCount
    fscanf(paramsFile, "%s %d", temp, &params->runsCount);
    //unused
    fscanf(paramsFile, "%s", temp);

    return params;
}


FdtdField  *initHostFieldWithParams(FdtdParams *)
{
    return NULL;
}


void deallocHostField(FdtdField *field)
{
}


FdtdField  *initDeviceFieldWithParams(FdtdParams *)
{
    return NULL;
}


void deallocDeviceField(FdtdField *field)
{
}
