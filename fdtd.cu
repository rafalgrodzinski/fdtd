#include "fdtd.h"

#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

#include "utils.h"


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
    deviceField = initDeviceFieldWithParams(params);

    printf("Reading materials data...\n");
    loadMaterials(params, field, "data/mat_specs_riken", params->inputPath);
    printf("Setting up mur boundary...\n");
    setupMurBoundary(params, field);

    // Setup CUDA parameters
    dim3 gridSize = dim3((params->nx + BLOCK_X - 1)/BLOCK_X,
                         (params->ny + BLOCK_Y - 1)/BLOCK_Y,
                         (params->nz + BLOCK_Z - 1)/BLOCK_Z);
    dim3 blockSize = dim3(BLOCK_X, BLOCK_Y, BLOCK_Z);

    // Main loop
    for(int i=0; i<params->iterationsCount; i++) {
        // Run 0
        printf("Running iteration %d", i);

        // Run 1
        printf("Running iteration %d", i+1);

        // Run 2
        printf("Running iteration %d", i+2);
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
    fscanf(paramsFile, "%s %f\n", temp, &params->eps_s);
    //feps_inf (eps_i)
    fscanf(paramsFile, "%s %f\n", temp, &params->eps_i);
    //ftau_d (tau_d)
    fscanf(paramsFile, "%s %f\n", temp, &params->tau_d);
    
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
    printf("Default eps_s:              %g\n", params->eps_s);
    printf("Default eps_i:              %g\n", params->eps_i);
    printf("Default tau_d:              %g\n", params->tau_d);
    printf("\n");
}


FdtdField *initFieldWithParams(FdtdParams *params)
{
    int n = params->nx * params->ny * params->nz * sizeof(float); 

    FdtdField *field = (FdtdField *)malloc(sizeof(FdtdField));

    // e
    CHECK(cudaHostAlloc(&field->ex1, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ey1, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ez1, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->ex2, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ey2, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ez2, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->ex3, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ey3, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ez3, n, cudaHostAllocDefault))

    // h
    CHECK(cudaHostAlloc(&field->hx, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->hy, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->hz, n, cudaHostAllocDefault))

    // d
    CHECK(cudaHostAlloc(&field->dx1, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dy1, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dz1, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->dx2, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dy2, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dz2, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->dx3, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dy3, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dz3, n, cudaHostAllocDefault))

    // sigma, eps, tau
    CHECK(cudaHostAlloc(&field->sigma, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->eps_s, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->eps_i, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->tau_d, n, cudaHostAllocDefault))

    // rp
    CHECK(cudaHostAlloc(&field->rp_x_0, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->rp_y_0, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->rp_z_0, n, cudaHostAllocDefault))

    CHECK(cudaHostAlloc(&field->rp_x_end, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->rp_y_end, n, cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->rp_z_end, n, cudaHostAllocDefault))

    return field;
}


void deallocField(FdtdField *field)
{
    // e
    CHECK(cudaFree(field->ex1))
    CHECK(cudaFree(field->ey1))
    CHECK(cudaFree(field->ez1))

    CHECK(cudaFree(field->ex2))
    CHECK(cudaFree(field->ey2))
    CHECK(cudaFree(field->ez2))

    CHECK(cudaFree(field->ex3))
    CHECK(cudaFree(field->ey3))
    CHECK(cudaFree(field->ez3))

    // h
    CHECK(cudaFree(field->hx))
    CHECK(cudaFree(field->hy))
    CHECK(cudaFree(field->hz))

    // d
    CHECK(cudaFree(field->dx1))
    CHECK(cudaFree(field->dy1))
    CHECK(cudaFree(field->dz1))

    CHECK(cudaFree(field->dx2))
    CHECK(cudaFree(field->dy2))
    CHECK(cudaFree(field->dz2))

    CHECK(cudaFree(field->dx3))
    CHECK(cudaFree(field->dy3))
    CHECK(cudaFree(field->dz3))

    // sigma, eps, tau
    CHECK(cudaFree(&field->sigma))
    CHECK(cudaFree(&field->eps_s))
    CHECK(cudaFree(&field->eps_i))
    CHECK(cudaFree(&field->tau_d))

    // rp
    CHECK(cudaFree(&field->rp_x_0))
    CHECK(cudaFree(&field->rp_y_0))
    CHECK(cudaFree(&field->rp_z_0))

    CHECK(cudaFree(&field->rp_x_end))
    CHECK(cudaFree(&field->rp_y_end))
    CHECK(cudaFree(&field->rp_z_end))

    free(field);
}


FdtdField *initDeviceFieldWithParams(FdtdParams *params)
{
    int n = params->nx * params->ny * params->nz * sizeof(float); 

    FdtdField *field = (FdtdField *)malloc(sizeof(FdtdField));

    // e
    CHECK(cudaMalloc(&field->ex1, n))
    CHECK(cudaMalloc(&field->ey1, n))
    CHECK(cudaMalloc(&field->ez1, n))

    CHECK(cudaMalloc(&field->ex2, n))
    CHECK(cudaMalloc(&field->ey2, n))
    CHECK(cudaMalloc(&field->ez2, n))

    CHECK(cudaMalloc(&field->ex3, n))
    CHECK(cudaMalloc(&field->ey3, n))
    CHECK(cudaMalloc(&field->ez3, n))

    // h
    CHECK(cudaMalloc(&field->hx, n))
    CHECK(cudaMalloc(&field->hy, n))
    CHECK(cudaMalloc(&field->hz, n))

    // d
    CHECK(cudaMalloc(&field->dx1, n))
    CHECK(cudaMalloc(&field->dy1, n))
    CHECK(cudaMalloc(&field->dz1, n))

    CHECK(cudaMalloc(&field->dx2, n))
    CHECK(cudaMalloc(&field->dy2, n))
    CHECK(cudaMalloc(&field->dz2, n))

    CHECK(cudaMalloc(&field->dx3, n))
    CHECK(cudaMalloc(&field->dy3, n))
    CHECK(cudaMalloc(&field->dz3, n))

    // sigma, eps, tau
    CHECK(cudaMalloc(&field->eps_i, n))
    CHECK(cudaMalloc(&field->eps_s, n))
    CHECK(cudaMalloc(&field->tau_d, n))
    CHECK(cudaMalloc(&field->sigma, n))

    // rp
    CHECK(cudaMalloc(&field->rp_x_0, n))
    CHECK(cudaMalloc(&field->rp_y_0, n))
    CHECK(cudaMalloc(&field->rp_z_0, n))

    CHECK(cudaMalloc(&field->rp_x_end, n))
    CHECK(cudaMalloc(&field->rp_y_end, n))
    CHECK(cudaMalloc(&field->rp_z_end, n))

    return field;
}


void deallocDeviceField(FdtdField *field)
{
    // e
    CHECK(cudaFree(field->ex1))
    CHECK(cudaFree(field->ey1))
    CHECK(cudaFree(field->ez1))

    CHECK(cudaFree(field->ex2))
    CHECK(cudaFree(field->ey2))
    CHECK(cudaFree(field->ez2))

    CHECK(cudaFree(field->ex3))
    CHECK(cudaFree(field->ey3))
    CHECK(cudaFree(field->ez3))

    // h
    CHECK(cudaFree(field->hx))
    CHECK(cudaFree(field->hy))
    CHECK(cudaFree(field->hz))

    // d
    CHECK(cudaFree(field->dx1))
    CHECK(cudaFree(field->dy1))
    CHECK(cudaFree(field->dz1))

    CHECK(cudaFree(field->dx2))
    CHECK(cudaFree(field->dy2))
    CHECK(cudaFree(field->dz2))

    CHECK(cudaFree(field->dx3))
    CHECK(cudaFree(field->dy3))
    CHECK(cudaFree(field->dz3))

    // sigma, eps, tau
    CHECK(cudaFree(field->eps_i))
    CHECK(cudaFree(field->eps_s))
    CHECK(cudaFree(field->tau_d))
    CHECK(cudaFree(field->sigma))

    // rp
    CHECK(cudaFree(field->rp_x_0))
    CHECK(cudaFree(field->rp_y_0))
    CHECK(cudaFree(field->rp_z_0))

    CHECK(cudaFree(field->rp_x_end))
    CHECK(cudaFree(field->rp_y_end))
    CHECK(cudaFree(field->rp_z_end))
}


void loadMaterials(FdtdParams *params, FdtdField *field, const char *specsFilePath, const char *materialsPath)
{
    // Load material specs
    int specsCount = 94;
    float *specs = (float *)calloc(specsCount * 4, sizeof(float));
    char temp[1024];
    int index;
    float sigma_value, eps_s_value, eps_i_value, tau_d_value;

    FILE *specsFile = fopen(specsFilePath, "r");
    //check(specsFile != NULL, "Cannot open specs file");

    for(int i=0; i<94; i++) {
        fscanf(specsFile, "%d %s %g %g %g %g", &index, temp, &sigma_value, &eps_s_value, &eps_i_value, &tau_d_value);
        specs[index*4 + 0] = sigma_value;
        specs[index*4 + 1] = eps_s_value;
        specs[index*4 + 2] = eps_i_value;
        specs[index*4 + 3] = tau_d_value;
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
                field->eps_s[offset] = specs[code*4 + 1];
                field->eps_i[offset] = specs[code*4 + 2];
                field->tau_d[offset] = specs[code*4 + 3];
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

                field->rp_x_0[offset] = 0.0; 
            }

            for(int ix = params->nx - 2; ix < params->nx; ix++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rp_x_end[offset] = 0.0;
            }
        }
    }

    // Setup rp_y
    for(int iz = 0; iz < params->nz; iz++) {
        for(int ix = 0; ix < params->nx; ix++) {
            for(int iy = 0; iy < 2; iy++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rp_y_0[offset] = 0.0; 
            }

            for(int iy = params->ny - 2; iy < params->ny; iy++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rp_y_end[offset] = 0.0;
            }
        }
    }

    // Setup rp_z
    for(int iy = 0; iy < params->ny; iy++) {
        for(int ix = 0; ix < params->nx; ix++) {
            for(int iz = 0; iz < 2; iz++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rp_z_0[offset] = 0.0; 
            }

            for(int iz = params->nz - 2; iz < params->nz; iz++) {
                int offset = iz * params->nx * params->ny + iy * params->nx + ix;

                field->rp_z_end[offset] = 0.0;
            }
        }
    }
}
