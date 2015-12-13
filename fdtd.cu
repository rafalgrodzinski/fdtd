#include "fdtd.h"

#include <stdlib.h>
#include <stdio.h>

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
    char temp[1024];
    int index;
    float sigmaValue, epsSValue, epsIValue, tauDValue;

    FILE *specsFile = fopen(specsFilePath, "r");
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
    // Setup rp_x
    for(int iz = 0; iz < params->nz; iz++) {
        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix = 0; ix < 2; ix++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rpx0[offset] = 0.0; 
            }

            for(int ix = params->nx - 2; ix < params->nx; ix++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rpxEnd[offset] = 0.0;
            }
        }
    }

    // Setup rp_y
    for(int iz = 0; iz < params->nz; iz++) {
        for(int ix = 0; ix < params->nx; ix++) {
            for(int iy = 0; iy < 2; iy++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rpy0[offset] = 0.0; 
            }

            for(int iy = params->ny - 2; iy < params->ny; iy++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rpyEnd[offset] = 0.0;
            }
        }
    }

    // Setup rp_z
    for(int iy = 0; iy < params->ny; iy++) {
        for(int ix = 0; ix < params->nx; ix++) {
            for(int iz = 0; iz < 2; iz++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rpz0[offset] = 0.0; 
            }

            for(int iz = params->nz - 2; iz < params->nz; iz++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rpzEnd[offset] = 0.0;
            }
        }
    }
}


void setupSources(FdtdParams *params, FdtdField *field)
{
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
