#include "fdtd.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "utils.h"
#include "fdtd_calculations.h"


#define BLOCK_X 4
#define BLOCK_Y 4
#define BLOCK_Z 4


int main(int argc, char **argv)
{
    // Read params
    FdtdParams *params;
    printf("Reading parameters...\n");
    params = initParamsWithPath("data/input_params");
    printParams(params);

    // Initialize field
    FdtdField  *field, *deviceField; // Used for CUDA

    printf("Initializing field...\n");
    field = initFieldWithParams(params);
    setupMurBoundary(params, field);
    printf("Initializing device field...\n");
    deviceField = initDeviceFieldWithParams(params);

    printf("Reading materials data...\n");
    loadMaterials(params, field, "data/mat_specs_riken", params->inputPath);

    printf("Initializing sources...\n");
    setupSources(params, field);

    printf("Copying data to GPU...\n");
    copyData(params, field, deviceField);

    // Setup CUDA parameters
    dim3 gridSize = dim3((params->nx + BLOCK_X - 1)/BLOCK_X,
                         (params->ny + BLOCK_Y - 1)/BLOCK_Y,
                         (params->nz + BLOCK_Z - 1)/BLOCK_Z);
    dim3 blockSize = dim3(BLOCK_X, BLOCK_Y, BLOCK_Z);

    // Main loop
    for(int i=0; i<params->iterationsCount; i += 3) {
        // Run 0
        printf("Running iteration %d", i);

        updateHField<<<gridSize, blockSize>>>(deviceField->hx,  deviceField->hy,  deviceField->hz,                    
                                              deviceField->ex2, deviceField->ey2, deviceField->ez2,                 
                                              params->nx, params->ny, params->nz,                 
                                              params->dt, params->dx, params->dy, params->dz, 
                                              params->mu0);
        CHECK(cudaDeviceSynchronize())

        updateDField<<<gridSize, blockSize>>>(deviceField->dx0, deviceField->dy0, deviceField->dz0, 
                                              deviceField->dx2, deviceField->dy2, deviceField->dz2, 
                                              deviceField->hx,  deviceField->hy,  deviceField->hz,    
                                              params->nx, params->ny, params->nz, 
                                              params->dt, params->dx, params->dy, params->dz);
        CHECK(cudaDeviceSynchronize())
 
        updateSources<<<gridSize, blockSize>>>(deviceField->dz0, deviceField->dz2,                                 
                                               deviceField->hx,  deviceField->hy,                                   
                                               params->nx, params->ny, params->nz,
                                               params->dt, params->dx, params->dy, params->dz, 
                                               params->sources, params->jz,                                
                                               params->sourcesCount, i);
        CHECK(cudaDeviceSynchronize())
            
        updateEField<<<gridSize, blockSize>>>(deviceField->ex0, deviceField->ey0, deviceField->ez0, 
                                              deviceField->ex2, deviceField->ey2, deviceField->ez2, 
                                              deviceField->ex1, deviceField->ey1, deviceField->ez1, 
                                              deviceField->dx0, deviceField->dy0, deviceField->dz0, 
                                              deviceField->dx2, deviceField->dy2, deviceField->dz2, 
                                              deviceField->dx1, deviceField->dy1, deviceField->dz1, 
                                              deviceField->sigma, deviceField->epsI, deviceField->epsS, deviceField->tauD,             
                                              params->nx, params->ny, params->nz, 
                                              params->dt, params->eps0);
        CHECK(cudaDeviceSynchronize())
            
        updateMurBoundary<<<gridSize, blockSize>>>(deviceField->ex0,  deviceField->ey0,  deviceField->ez0,                 
                                                   deviceField->ex2,  deviceField->ey2,  deviceField->ez2,                 
                                                   deviceField->rpx0, deviceField->rpy0, deviceField->rpz0,                         
                                                   deviceField->rpxEnd, deviceField->rpyEnd, deviceField->rpzEnd,                         
                                                   params->nx, params->ny, params->nz,                 
                                                   params->dt, params->dx, params->dy, params->dz, 
                                                   params->mu0, params->eps0);
        CHECK(cudaDeviceSynchronize());

        // Write results
        writeResults(params, field,
                     field->ex0, field->ey0, field->ez0,
                     field->dx0, field->dy0, field->dz0,
                     i, params->outputPath);

        // Run 1
        printf("Running iteration %d", i+1);

        updateHField<<<gridSize, blockSize>>>(deviceField->hx,  deviceField->hy,  deviceField->hz,                    
                                              deviceField->ex0, deviceField->ey0, deviceField->ez0,                 
                                              params->nx, params->ny, params->nz,                 
                                              params->dt, params->dx, params->dy, params->dz, 
                                              params->mu0);
        CHECK(cudaDeviceSynchronize())

        updateDField<<<gridSize, blockSize>>>(deviceField->dx1, deviceField->dy1, deviceField->dz1, 
                                              deviceField->dx0, deviceField->dy0, deviceField->dz0, 
                                              deviceField->hx,  deviceField->hy,  deviceField->hz,    
                                              params->nx, params->ny, params->nz, 
                                              params->dt, params->dx, params->dy, params->dz);
        CHECK(cudaDeviceSynchronize())
 
        updateSources<<<gridSize, blockSize>>>(deviceField->dz1, deviceField->dz0,                                 
                                               deviceField->hx,  deviceField->hy,                                   
                                               params->nx, params->ny, params->nz,
                                               params->dt, params->dx, params->dy, params->dz, 
                                               params->sources, params->jz,                                
                                               params->sourcesCount, i);
        CHECK(cudaDeviceSynchronize())
            
        updateEField<<<gridSize, blockSize>>>(deviceField->ex1, deviceField->ey1, deviceField->ez1, 
                                              deviceField->ex0, deviceField->ey0, deviceField->ez0, 
                                              deviceField->ex2, deviceField->ey2, deviceField->ez2, 
                                              deviceField->dx1, deviceField->dy1, deviceField->dz1, 
                                              deviceField->dx0, deviceField->dy0, deviceField->dz0, 
                                              deviceField->dx2, deviceField->dy2, deviceField->dz2, 
                                              deviceField->sigma, deviceField->epsI, deviceField->epsS, deviceField->tauD,             
                                              params->nx, params->ny, params->nz, 
                                              params->dt, params->eps0);
        CHECK(cudaDeviceSynchronize())
            
        updateMurBoundary<<<gridSize, blockSize>>>(deviceField->ex1,  deviceField->ey1,  deviceField->ez1,                 
                                                   deviceField->ex0,  deviceField->ey0,  deviceField->ez0,                 
                                                   deviceField->rpx0, deviceField->rpy0, deviceField->rpz0,                         
                                                   deviceField->rpxEnd, deviceField->rpyEnd, deviceField->rpzEnd,                         
                                                   params->nx, params->ny, params->nz,                 
                                                   params->dt, params->dx, params->dy, params->dz, 
                                                   params->mu0, params->eps0);
        CHECK(cudaDeviceSynchronize())

        // Write results
        writeResults(params, field,
                     field->ex1, field->ey1, field->ez1,
                     field->dx1, field->dy1, field->dz1,
                     i+1, params->outputPath);

        // Run 2
        printf("Running iteration %d", i+2);

        updateHField<<<gridSize, blockSize>>>(deviceField->hx,  deviceField->hy,  deviceField->hz,                    
                                              deviceField->ex1, deviceField->ey1, deviceField->ez1,                 
                                              params->nx, params->ny, params->nz,                 
                                              params->dt, params->dx, params->dy, params->dz, 
                                              params->mu0);
        CHECK(cudaDeviceSynchronize())

        updateDField<<<gridSize, blockSize>>>(deviceField->dx2, deviceField->dy2, deviceField->dz2, 
                                              deviceField->dx1, deviceField->dy1, deviceField->dz1, 
                                              deviceField->hx,  deviceField->hy,  deviceField->hz,    
                                              params->nx, params->ny, params->nz, 
                                              params->dt, params->dx, params->dy, params->dz);
        CHECK(cudaDeviceSynchronize())
 
        updateSources<<<gridSize, blockSize>>>(deviceField->dz2, deviceField->dz1,                                 
                                               deviceField->hx,  deviceField->hy,                                   
                                               params->nx, params->ny, params->nz,
                                               params->dt, params->dx, params->dy, params->dz, 
                                               params->sources, params->jz,                                
                                               params->sourcesCount, i);
        CHECK(cudaDeviceSynchronize())
            
        updateEField<<<gridSize, blockSize>>>(deviceField->ex2, deviceField->ey2, deviceField->ez2, 
                                              deviceField->ex1, deviceField->ey1, deviceField->ez1, 
                                              deviceField->ex0, deviceField->ey0, deviceField->ez0, 
                                              deviceField->dx2, deviceField->dy2, deviceField->dz2, 
                                              deviceField->dx1, deviceField->dy1, deviceField->dz1, 
                                              deviceField->dx0, deviceField->dy0, deviceField->dz0, 
                                              deviceField->sigma, deviceField->epsI, deviceField->epsS, deviceField->tauD,             
                                              params->nx, params->ny, params->nz, 
                                              params->dt, params->eps0);
        CHECK(cudaDeviceSynchronize())
            
        updateMurBoundary<<<gridSize, blockSize>>>(deviceField->ex2,  deviceField->ey2,  deviceField->ez2,                 
                                                   deviceField->ex1,  deviceField->ey1,  deviceField->ez1,                 
                                                   deviceField->rpx0, deviceField->rpy0, deviceField->rpz0,                         
                                                   deviceField->rpxEnd, deviceField->rpyEnd, deviceField->rpzEnd,                         
                                                   params->nx, params->ny, params->nz,                 
                                                   params->dt, params->dx, params->dy, params->dz, 
                                                   params->mu0, params->eps0);
        CHECK(cudaDeviceSynchronize())

        // Write results
        writeResults(params, field,
                     field->ex2, field->ey2, field->ez2,
                     field->dx2, field->dy2, field->dz2,
                     i+2, params->outputPath);
    }

    // Clean up
    deallocDeviceField(deviceField);
    deallocField(field);
    deallocParams(params);
}


FdtdParams *initParamsWithPath(const char *filePath)
{
    FdtdParams *params = (FdtdParams *)malloc(sizeof(FdtdParams));
    params->inputPath = (char *)malloc(sizeof(char) * 1024);
    params->outputPath = (char *)malloc(sizeof(char) * 1024);

    FILE *paramsFile = fopen(filePath, "r");
    //check(paramsFile != NULL, "Cannot open file");
    
    int tempLength = 1024;
    char temp[tempLength];

    //nx ny nz (field size)
    fscanf(paramsFile, "%s %d %d %d\n", temp, &params->nx, &params->ny, &params->nz);
    //t_max (simulation runs count)
    fscanf(paramsFile, "%s %d\n", temp, &params->iterationsCount);
    params->iterationsCount = ((params->iterationsCount - 1)/3 + 1) * 3; // Has to be divisible by 3
    //unused (nf)
    fgets(temp, tempLength, paramsFile);
    //env_set_dir (input path)
    fscanf(paramsFile, "%s %s\n", temp, params->inputPath);
    //unused (env_file_prefix)
    fgets(temp, tempLength, paramsFile);
    //output_dir (output path) 
    fscanf(paramsFile, "%s %s\n", temp, params->outputPath);
    //unused (output_format)
    fgets(temp, tempLength, paramsFile);
    //unused (impulse_resp_flag)
    fgets(temp, tempLength, paramsFile);
    //unused (pec_flag) 
    fgets(temp, tempLength, paramsFile);
    //unused (read_env_flag)
    fgets(temp, tempLength, paramsFile);
    //unused (output_flag)
    fgets(temp, tempLength, paramsFile);
    //unused (bzip2_flag)
    fgets(temp, tempLength, paramsFile);
    //unused (output_start)
    fgets(temp, tempLength, paramsFile);
    //unused (output_finish)
    fgets(temp, tempLength, paramsFile);
    //unused (source_type)
    fgets(temp, tempLength, paramsFile);
    //elements_per_wavelength
    fscanf(paramsFile, "%s %d\n", temp, &params->elementsPerWave);
    //wave_freq
    fscanf(paramsFile, "%s %g\n", temp, &params->waveFrequency);
    //pulse_width
    fscanf(paramsFile, "%s %g\n", temp, &params->pulseWidth);
    //pulse_modulation_frequency
    fscanf(paramsFile, "%s %g\n", temp, &params->pulseModulationFrequency);
    //number_of_excitation_sources
    fscanf(paramsFile, "%s %d\n", temp, &params->sourcesCount);
    //source_location
    params->sources = (int *)malloc(sizeof(int) * params->sourcesCount * 3);
    for(int i=0; i<params->sourcesCount; i++) {
        fscanf(paramsFile, "%s %d %d %d\n", temp,
                                            &params->sources[i*3],
                                            &params->sources[i*3 + 1],
                                            &params->sources[i*3 + 2]);
    }
    //unused (pulse_type)
    fgets(temp, tempLength, paramsFile);
    //fsigma (sigma)
    fscanf(paramsFile, "%s %f\n", temp, &params->sigma);
    //feps_s (eps_s)
    fscanf(paramsFile, "%s %f\n", temp, &params->epsS);
    //feps_inf (eps_i)
    fscanf(paramsFile, "%s %f\n", temp, &params->epsI);
    //ftau_d (tau_d)
    fscanf(paramsFile, "%s %f\n", temp, &params->tauD);
    
    fclose(paramsFile);

    // Generate rest of the values
    params->pi = acos(-1.0);
    params->c = 3.0 * pow(10.0, 8.0);
    params->timeskip = 1.0;
    params->lambda = params->c / params->waveFrequency;
    params->dx = params->lambda / params->elementsPerWave;
    params->dy = params->dx;
    params->dz = params->dx;
    params->dt = 1.0 * params->timeskip / (params->c * sqrt(1.0/pow(params->dx, 2.0) + 1.0/pow(params->dy, 2.0) + 1.0/pow(params->dz, 2.0)));
    params->mu0 = 4.0 * params->pi * pow(10.0, -7.0);
    params->eps0 = 1.0 / params->mu0 * (params->c * params->c);

    return params;
}


void deallocParams(FdtdParams *params)
{
    free(params->inputPath);
    free(params->outputPath);
    free(params);
}


void printParams(FdtdParams *params)
{
    printf("Field size:                 %dx%dx%d\n", params->nx, params->ny, params->nz);
    printf("Iterations count:           %d\n", params->iterationsCount);
    printf("Input path:                 %s\n", params->inputPath);
    printf("Output path:                %s\n", params->outputPath);
    printf("Elements per wavelength:    %d\n", params->elementsPerWave);
    printf("Wave frequency:             %g\n", params->waveFrequency);
    printf("Pulse width:                %g\n", params->pulseWidth);
    printf("Pulse modulation frequency: %g\n", params->pulseModulationFrequency);
    printf("Sources count:              %d\n", params->sourcesCount);
    for(int i=0; i<params->sourcesCount; i++)
        printf("Source position:            %dx%dx%d\n", params->sources[i*3],
                                                         params->sources[i*3 + 1],
                                                         params->sources[i*3 + 2]);
    printf("Default sigma:              %g\n", params->sigma);
    printf("Default eps_s:              %g\n", params->epsS);
    printf("Default eps_i:              %g\n", params->epsI);
    printf("Default tau_d:              %g\n", params->tauD);
    printf("\n");
}


FdtdField *initFieldWithParams(FdtdParams *params)
{
    int n = params->nx * params->ny * params->nz * sizeof(float); 

    FdtdField *field = (FdtdField *)malloc(sizeof(FdtdField));

    // e
    CHECK(cudaHostAlloc(&field->ex0, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ey0, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ez0, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->ex1, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ey1, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ez1, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->ex2, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ey2, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ez2, n, cudaHostAllocDefault))

    // h
    CHECK(cudaHostAlloc(&field->hx, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->hy, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->hz, n, cudaHostAllocDefault))

    // d
    CHECK(cudaHostAlloc(&field->dx0, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dy0, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dz0, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->dx1, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dy1, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dz1, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->dx2, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dy2, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dz2, n, cudaHostAllocDefault))

    // sigma, eps, tau
    CHECK(cudaHostAlloc(&field->sigma, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->epsS, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->epsI, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->tauD, n, cudaHostAllocDefault))

    // rp
    CHECK(cudaHostAlloc(&field->rpx0, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->rpy0, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->rpz0, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->rpxEnd, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->rpyEnd, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->rpzEnd, n, cudaHostAllocDefault))

    return field;
}


void deallocField(FdtdField *field)
{
    // e
    CHECK(cudaFree(field->ex0))
    CHECK(cudaFree(field->ey0))
    CHECK(cudaFree(field->ez0))

    CHECK(cudaFree(field->ex1))
    CHECK(cudaFree(field->ey1))
    CHECK(cudaFree(field->ez1))

    CHECK(cudaFree(field->ex2))
    CHECK(cudaFree(field->ey2))
    CHECK(cudaFree(field->ez2))

    // h
    CHECK(cudaFree(field->hx))
    CHECK(cudaFree(field->hy))
    CHECK(cudaFree(field->hz))

    // d
    CHECK(cudaFree(field->dx0))
    CHECK(cudaFree(field->dy0))
    CHECK(cudaFree(field->dz0))

    CHECK(cudaFree(field->dx1))
    CHECK(cudaFree(field->dy1))
    CHECK(cudaFree(field->dz1))

    CHECK(cudaFree(field->dx2))
    CHECK(cudaFree(field->dy2))
    CHECK(cudaFree(field->dz2))

    // sigma, eps, tau
    CHECK(cudaFree(&field->sigma))
    CHECK(cudaFree(&field->epsS))
    CHECK(cudaFree(&field->epsI))
    CHECK(cudaFree(&field->tauD))

    // rp
    CHECK(cudaFree(&field->rpx0))
    CHECK(cudaFree(&field->rpy0))
    CHECK(cudaFree(&field->rpz0))

    CHECK(cudaFree(&field->rpxEnd))
    CHECK(cudaFree(&field->rpyEnd))
    CHECK(cudaFree(&field->rpzEnd))

    free(field);
}


FdtdField *initDeviceFieldWithParams(FdtdParams *params)
{
    int n = params->nx * params->ny * params->nz * sizeof(float); 

    FdtdField *field = (FdtdField *)malloc(sizeof(FdtdField));

    // e
    CHECK(cudaMalloc(&field->ex0, n))
    CHECK(cudaMalloc(&field->ey0, n))
    CHECK(cudaMalloc(&field->ez0, n))

    CHECK(cudaMalloc(&field->ex1, n))
    CHECK(cudaMalloc(&field->ey1, n))
    CHECK(cudaMalloc(&field->ez1, n))

    CHECK(cudaMalloc(&field->ex2, n))
    CHECK(cudaMalloc(&field->ey2, n))
    CHECK(cudaMalloc(&field->ez2, n))

    // h
    CHECK(cudaMalloc(&field->hx, n))
    CHECK(cudaMalloc(&field->hy, n))
    CHECK(cudaMalloc(&field->hz, n))

    // d
    CHECK(cudaMalloc(&field->dx0, n))
    CHECK(cudaMalloc(&field->dy0, n))
    CHECK(cudaMalloc(&field->dz0, n))

    CHECK(cudaMalloc(&field->dx1, n))
    CHECK(cudaMalloc(&field->dy1, n))
    CHECK(cudaMalloc(&field->dz1, n))

    CHECK(cudaMalloc(&field->dx2, n))
    CHECK(cudaMalloc(&field->dy2, n))
    CHECK(cudaMalloc(&field->dz2, n))

    // sigma, eps, tau
    CHECK(cudaMalloc(&field->epsI, n))
    CHECK(cudaMalloc(&field->epsS, n))
    CHECK(cudaMalloc(&field->tauD, n))
    CHECK(cudaMalloc(&field->sigma, n))

    // rp
    CHECK(cudaMalloc(&field->rpx0, n))
    CHECK(cudaMalloc(&field->rpy0, n))
    CHECK(cudaMalloc(&field->rpz0, n))

    CHECK(cudaMalloc(&field->rpxEnd, n))
    CHECK(cudaMalloc(&field->rpyEnd, n))
    CHECK(cudaMalloc(&field->rpzEnd, n))

    return field;
}


void deallocDeviceField(FdtdField *field)
{
    // e
    CHECK(cudaFree(field->ex0))
    CHECK(cudaFree(field->ey0))
    CHECK(cudaFree(field->ez0))

    CHECK(cudaFree(field->ex1))
    CHECK(cudaFree(field->ey1))
    CHECK(cudaFree(field->ez1))

    CHECK(cudaFree(field->ex2))
    CHECK(cudaFree(field->ey2))
    CHECK(cudaFree(field->ez2))

    // h
    CHECK(cudaFree(field->hx))
    CHECK(cudaFree(field->hy))
    CHECK(cudaFree(field->hz))

    // d
    CHECK(cudaFree(field->dx0))
    CHECK(cudaFree(field->dy0))
    CHECK(cudaFree(field->dz0))

    CHECK(cudaFree(field->dx1))
    CHECK(cudaFree(field->dy1))
    CHECK(cudaFree(field->dz1))

    CHECK(cudaFree(field->dx2))
    CHECK(cudaFree(field->dy2))
    CHECK(cudaFree(field->dz2))

    // sigma, eps, tau
    CHECK(cudaFree(field->epsI))
    CHECK(cudaFree(field->epsS))
    CHECK(cudaFree(field->tauD))
    CHECK(cudaFree(field->sigma))

    // rp
    CHECK(cudaFree(field->rpx0))
    CHECK(cudaFree(field->rpy0))
    CHECK(cudaFree(field->rpz0))

    CHECK(cudaFree(field->rpxEnd))
    CHECK(cudaFree(field->rpyEnd))
    CHECK(cudaFree(field->rpzEnd))
}


void loadMaterials(FdtdParams *params, FdtdField *field, const char *specsFilePath, const char *materialsPath)
{
    // Load material specs
    int specsCount = 94;
    float *specs = (float *)calloc(specsCount * 4, sizeof(float));
    if(specs == NULL) {
        printf("Couldn't alocate %d bytes for specs\n", specsCount*4*sizeof(float));
        exit(EXIT_FAILURE);
    }
    char temp[1024];
    int index;
    float sigmaValue, epsSValue, epsIValue, tauDValue;

    FILE *specsFile = fopen(specsFilePath, "r");
    if(specsFile == NULL) {
        printf("Couldn\'t open file %s\n", specsFilePath);
        exit(EXIT_FAILURE);
    }
    //check(specsFile != NULL, "Cannot open specs file");

    for(int i=0; i<94; i++) {
        fscanf(specsFile, "%d %s %g %g %g %g", &index, temp, &sigmaValue, &epsSValue, &epsIValue, &tauDValue);
        specs[index*4 + 0] = sigmaValue;
        specs[index*4 + 1] = epsSValue;
        specs[index*4 + 2] = epsIValue;
        specs[index*4 + 3] = tauDValue;
    }

    fclose(specsFile);

    // Load materials
    for(int iz=0; iz<params->nz; iz++) {
        char materialFileName[1024];
        sprintf(materialFileName, "%s/v1_%5d.pgm", materialsPath, iz);
        FILE *materialFile = fopen(materialFileName, "r");

        if(materialFile == NULL) {
            printf("Couldn\'t open file %s\n", materialFileName);
            exit(EXIT_FAILURE);
        }

        printf("Reading %s...", materialFileName);

        int width, height;
        fscanf(materialFile, "%s %s %s %d %d %s", temp, temp, temp, &width, &height, temp);

        for(int iy=0; iy<params->ny; iy++) {
            for(int ix=0; ix<params->nx; ix++) {
                int code;
                fscanf(materialFile, "%d", &code);

                int offset = iz*params->nx*params->ny + iy*params->nx + ix;
                field->sigma[offset] = specs[code*4 + 0];
                field->epsS[offset] = specs[code*4 + 1];
                field->epsI[offset] = specs[code*4 + 2];
                field->tauD[offset] = specs[code*4 + 3];
            }
        }

        fclose(materialFile);
    }

    free(specs);
}


void setupMurBoundary(FdtdParams *params, FdtdField *field)
{
#ifndef __APPLE__
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    // Setup rpx
    for(int iz = 0; iz < nz; iz++) {
        for(int iy = 0; iy < ny; iy++) {
            for(int ix = 0; ix < 2; ix++) {
                float complex c1 = 0.0 + 2.0 * params->pi * params->waveFrequency * OFFSET(field->tauD, ix, iy,iz);
                float complex c2 = 0.0 + OFFSET(field->sigma, ix, iy, iz) / (2.0 * params->pi * params->waveFrequency * params->eps0);

                OFFSET(field->rpx0, ix, iy, iz) = creal(OFFSET(field->epsI, ix, iy, iz) +
                                                        (OFFSET(field->epsS, ix, iy, iz) - OFFSET(field->epsI, ix, iy, iz)) / (1.0 + c1) - c2);
            }

            for(int ix = nx - 2; ix < nx; ix++) {
                float complex c1 = 0.0 + 2.0 * params->pi * params->waveFrequency * OFFSET(field->tauD, ix, iy, iz);
                float complex c2 = 0.0 + OFFSET(field->sigma, ix, iy, iz) / (2.0 * params->pi * params->waveFrequency * params->eps0);
                
                OFFSET(field->rpxEnd, ix, iy, iz) = creal(OFFSET(field->epsI, ix, iy, iz) +                                                  
                                                          (OFFSET(field->epsS, ix, iy, iz) - OFFSET(field->epsI, ix, iy, iz)) / (1.0 + c1) - c2);
            }
        }
    }

    // Setup rpy
    for(int iz = 0; iz < nz; iz++) {
        for(int ix = 0; ix < nx; ix++) {
            for(int iy = 0; iy < 2; iy++) {
                float complex c1 = 0.0 + 2.0 * params->pi * params->waveFrequency * OFFSET(field->tauD, ix, iy, iz) * I;
                float complex c2 = 0.0 + OFFSET(field->sigma, ix, iy, iz) /(2.0 * params->pi * params->waveFrequency * params->eps0) * I;
                
                OFFSET(field->rpy0, ix, iy, iz) = creal(OFFSET(field->epsI, ix, iy, iz) +                                                      
                                                        (OFFSET(field->epsS, ix, iy, iz) - OFFSET(field->epsI, ix, iy, iz)) / (1.0 + c1) - c2);
            }

            for(int iy = ny - 2; iy < ny; iy++) {
                float complex c1 = 0.0 + 2.0 * params->pi * params->waveFrequency * OFFSET(field->tauD, ix, iy, iz) * I;
                float complex c2 = 0.0 + OFFSET(field->sigma, ix, iy, iz) / (2 * params->pi * params->waveFrequency * params->eps0) * I;
                
                OFFSET(field->rpyEnd, ix, iy, iz) = creal(OFFSET(field->epsI, ix, iy, iz) +                                                      
                                                          (OFFSET(field->epsS, ix, iy, iz) - OFFSET(field->epsI, ix, iy, iz)) / (1.0 + c1) - c2);
            }
        }
    }

    // Setup rpz
    for(int iy = 0; iy < ny; iy++) {
        for(int ix = 0; ix < nx; ix++) {
            for(int iz = 0; iz < 2; iz++) {
                float complex c1 = 0.0 + 2.0 * params->pi * params->waveFrequency * OFFSET(field->tauD, ix, iy, iz) * I;
                float complex c2 = 0.0 + OFFSET(field->sigma, ix, iy, iz) / (2.0 * params->pi * params->waveFrequency * params->eps0) * I;
                
                OFFSET(field->rpz0, ix, iy, iz) = creal(OFFSET(field->epsI, ix, iy, iz) +                                                  
                                                        (OFFSET(field->epsS, ix, iy, iz) - OFFSET(field->epsI, ix, iy, iz)) / (1.0 + c1) - c2);
            }

            for(int iz = nz - 2; iz < nz; iz++) {
                float complex c1 = 0.0 + 2.0 * params->pi * params->waveFrequency * OFFSET(field->tauD, ix, iy, iz) * I;
                float complex c2 = 0.0 + OFFSET(field->sigma, ix, iy, iz) / (2.0 * params->pi * params->waveFrequency * params->eps0) * I;
                
                OFFSET(field->rpzEnd, ix, iy, iz) = creal(OFFSET(field->epsI, ix, iy, iz) +                                                  
                                                          (OFFSET(field->epsS, ix, iy, iz) - OFFSET(field->epsI, ix, iy, iz)) / (1.0 + c1) - c2);
            }
        }
    }
#endif
}


void setupSources(FdtdParams *params, FdtdField *field)
{
    int fine;
    int temp;
    int i2;
    int istart;
    float *tmpdata, *tmpdata2;

    int tmpOff = 1<<16;
    params->jz = (float *)calloc(tmpOff, sizeof(float));
    tmpdata = (float *)calloc(tmpOff * 2, sizeof(float));
    tmpdata2 = (float *)calloc(tmpOff * 2, sizeof(float));
    
    // fine & temp
    fine = (1<<13) * params->pulseWidth * params->waveFrequency * params->dt;

    temp = 1.0/(params->pulseWidth * params->waveFrequency)/(params->dt / fine)/2.0;
    
    // tmpdata
    for(int i = -temp - 1; i <= temp + 1; i++) {
        float value = exp(pow(-((float)i/(((float)temp + 1.0)/4.0)), 2.0));
        value *= cos(2.0 * acos(-1.0) * params->pulseModulationFrequency * params->waveFrequency * i * (params->dt / fine));
        tmpdata[i + tmpOff] = value; 
    }

    // istart
    int off = 1<<12;
    for(int i = -off; i < off; i++) {
         if((fabs(tmpdata[i + tmpOff]) >= pow(1.0, -9.0)) && (i % fine == 0)) {
            istart = i;
            break;
         }
    }
    
    // setup jz 1/2
    i2 = 0;
    for(int i=istart; i <= temp+1; i++) {
          params->jz[i2] = tmpdata[i] * pow(10.0, -15.0) / params->dt / 3.0;
          i2++;
    }
    
    //setup tmpdata2
    for(int i=2; i <= (1<<14); i++) {
        tmpdata2[i - 1 + tmpOff] = (((params->jz[i + 1 - 1] - params->jz[i - 1]) / params->dt) +
                                    ((params->jz[i - 1] - params->jz[i - 1 - 1]) / params->dt)) / 
                                    2.0 * (params->dt * params->dz) / (params->dx * params->dy * params->dz);
    }
    
    // Setup jz 2/2
    for(int i=0; i < (2 << 14); i++) {
        params->jz[i] = tmpdata2[i + 1];
    }

    free(tmpdata2);
    free(tmpdata);
}


void copyData(FdtdParams *params, FdtdField *field, FdtdField *deviceField)
{
    int n = params->nx * params->ny * params->nz * sizeof(float); 

    //e
    cudaMemcpy(deviceField->ex0, field->ex0, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->ey0, field->ey0, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->ez0, field->ez0, n, cudaMemcpyHostToDevice);

    cudaMemcpy(deviceField->ex1, field->ex1, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->ey1, field->ey1, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->ez1, field->ez1, n, cudaMemcpyHostToDevice);

    cudaMemcpy(deviceField->ex2, field->ex2, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->ey2, field->ey2, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->ez2, field->ez2, n, cudaMemcpyHostToDevice);

    //h
    cudaMemcpy(deviceField->hx, field->hx, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->hy, field->hy, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->hz, field->hz, n, cudaMemcpyHostToDevice);

    //d
    cudaMemcpy(deviceField->dx0, field->dx0, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->dy0, field->dy0, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->dz0, field->dz0, n, cudaMemcpyHostToDevice);

    cudaMemcpy(deviceField->dx1, field->dx1, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->dy1, field->dy1, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->dz1, field->dz1, n, cudaMemcpyHostToDevice);

    cudaMemcpy(deviceField->dx2, field->dx2, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->dy2, field->dy2, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->dz2, field->dz2, n, cudaMemcpyHostToDevice);

    cudaMemcpy(deviceField->epsI, field->epsI, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->epsS, field->epsS, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->tauD, field->tauD, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->sigma, field->sigma, n, cudaMemcpyHostToDevice);

    cudaMemcpy(deviceField->rpx0, field->rpx0, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->rpy0, field->rpy0, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->rpz0, field->rpz0, n, cudaMemcpyHostToDevice);

    cudaMemcpy(deviceField->rpxEnd, field->rpxEnd, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->rpyEnd, field->rpyEnd, n, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceField->rpzEnd, field->rpzEnd, n, cudaMemcpyHostToDevice);
}


void writeResults(FdtdParams *params, FdtdField *field,
                  float *exSource, float *eySource, float *ezSource,
                  float *dxSource, float *dySource, float *dzSource,
                  int currentIteration, char *outputPath)
{
    char outputFilePath[1024];
    FILE *outputFile;

    // Used by OFFSET macro
    int nx = params->nx;
    int ny = params->ny;

    // Output x
    sprintf(outputPath, "%s/E_field_x_%5d.out", outputPath, currentIteration);

    outputFile = fopen(outputFilePath, "w");
    for(int isrc=0; isrc < params->sourcesCount; isrc++) {
        int iy = params->sources[isrc * 3 + 1];
        int iz = params->sources[isrc * 3 + 2];
        for(int ix=0; ix < params->nx; ix++) {
            fprintf(outputFile, "%4d %4d %4d %g %g %g %g %g %g %g %g %g", ix, iy, iz,
                    OFFSET(dxSource, ix, iy, iz),  OFFSET(dySource, ix, iy, iz),  OFFSET(dzSource, ix, iy, iz),
                    OFFSET(field->hx, ix, iy, iz), OFFSET(field->hy, ix, iy, iz), OFFSET(field->hz, ix, iy, iz),
                    OFFSET(exSource, ix, iy, iz),  OFFSET(eySource, ix, iy, iz),  OFFSET(ezSource, ix, iy, iz));
        }
    }
    fclose(outputFile);

    // Output y
    sprintf(outputPath, "%s/E_field_y_%5d.out", outputPath, currentIteration);

    outputFile = fopen(outputFilePath, "w");
    for(int isrc=0; isrc < params->sourcesCount; isrc++) {
        int ix = params->sources[isrc * 3 + 0];
        int iz = params->sources[isrc * 3 + 2];
        for(int iy=0; iy < params->ny; iy++) {
            fprintf(outputFile, "%4d %4d %4d %g %g %g %g %g %g %g %g %g", ix, iy, iz,
                    OFFSET(dxSource, ix, iy, iz),  OFFSET(dySource, ix, iy, iz),  OFFSET(dzSource, ix, iy, iz),
                    OFFSET(field->hx, ix, iy, iz), OFFSET(field->hy, ix, iy, iz), OFFSET(field->hz, ix, iy, iz),
                    OFFSET(exSource, ix, iy, iz),  OFFSET(eySource, ix, iy, iz),  OFFSET(ezSource, ix, iy, iz));
        }
    }
    fclose(outputFile);

    // Output z
    sprintf(outputPath, "%s/E_field_z_%5d.out", outputPath, currentIteration);

    outputFile = fopen(outputFilePath, "w");
    for(int isrc=0; isrc < params->sourcesCount; isrc++) {
        int ix = params->sources[isrc * 3 + 0];
        int iy = params->sources[isrc * 3 + 1];
        for(int iz=0; iz < params->nz; iz++) {
            fprintf(outputFile, "%4d %4d %4d %g %g %g %g %g %g %g %g %g", ix, iy, iz,
                    OFFSET(dxSource, ix, iy, iz),  OFFSET(dySource, ix, iy, iz),  OFFSET(dzSource, ix, iy, iz),
                    OFFSET(field->hx, ix, iy, iz), OFFSET(field->hy, ix, iy, iz), OFFSET(field->hz, ix, iy, iz),
                    OFFSET(exSource, ix, iy, iz),  OFFSET(eySource, ix, iy, iz),  OFFSET(ezSource, ix, iy, iz));
        }
    }
    fclose(outputFile);
}
