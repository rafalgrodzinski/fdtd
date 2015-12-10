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
    params = initParamsWithPath("data/input_params");
    printParams(params);

    // Initialize field
    FdtdField  *hostField, *deviceField; // Used for CUDA

    hostField = initHostFieldWithParams(params);
    deviceField = initDeviceFieldWithParams(params);

    // Setup CUDA parameters
    dim3 gridSize = dim3((params->nx + BLOCK_X - 1)/BLOCK_X,
                         (params->ny + BLOCK_Y - 1)/BLOCK_Y,
                         (params->nz + BLOCK_Z - 1)/BLOCK_Z);
    dim3 blockSize = dim3(BLOCK_X, BLOCK_Y, BLOCK_Z);

    // Main loop
    for(int i=0; i<params->iterationsCount; i++) {
        // Run 0
        // Run 1
        // Run 2
    }

    // Clean up
    deallocParams(params);
}


void printUsage()
{
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
}


FdtdField  *initHostFieldWithParams(FdtdParams *params)
{
    return NULL;
}


void deallocHostField(FdtdField *field)
{
}


FdtdField  *initDeviceFieldWithParams(FdtdParams *params)
{
    return NULL;
}


void deallocDeviceField(FdtdField *field)
{
}
