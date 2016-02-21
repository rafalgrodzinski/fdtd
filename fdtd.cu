#include "fdtd.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "utils.h"
#include "fdtd_calculations.h"


#define BLOCK_X 128
#define BLOCK_Y 1
#define BLOCK_Z 1


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
    setupSources(params);

    printf("Copying data to GPU...\n\n");
    copyDataToDevice(params, field, deviceField);
    copySymbolsToDevice(params);

    // Setup CUDA parameters
    dim3 gridSize = dim3((params->nx + BLOCK_X - 1)/BLOCK_X,
                         (params->ny + BLOCK_Y - 1)/BLOCK_Y,
                         (params->nz + BLOCK_Z - 1)/BLOCK_Z);
    dim3 blockSize = dim3(BLOCK_X, BLOCK_Y, BLOCK_Z);

    // Create streams
    cudaStream_t streamH;
    cudaStream_t streamD;
    cudaStream_t streamE;

    CHECK(cudaStreamCreate(&streamH))
    CHECK(cudaStreamCreate(&streamD))
    CHECK(cudaStreamCreate(&streamE))
    
    cudaEvent_t eventH;
    cudaEvent_t eventD;
    cudaEvent_t eventE;

    CHECK(cudaEventCreate(&eventH))
    CHECK(cudaEventCreate(&eventD))
    CHECK(cudaEventCreate(&eventE))
    
    int bytesCount = params->nx * params->ny * params->nz * sizeof(float); 

    // Threads
    CopyParams *hCopyParams;
    CopyParams *dCopyParams;
    CopyParams *eCopyParams;

    pthread_t *hThread = NULL;
    pthread_t *dThread = NULL;
    pthread_t *eThread = NULL;

    ResultsParams *resultsParams;
    pthread_t *threads = (pthread_t *)malloc(params->iterationsCount * sizeof(pthread_t));

    CHECK(cudaEventRecord(eventE))

    // Main loop
    for(int i=0; i<params->iterationsCount; i += 3) {

        // Run 0
        printf("Running iteration %d\n", i);

        // H field
        CHECK(cudaStreamWaitEvent(streamH, eventE, 0));

        if(hThread != NULL)
            pthread_join(*hThread, NULL);

        updateHField<<<gridSize, blockSize, 0, streamH>>>(deviceField->hx,  deviceField->hy,  deviceField->hz,                    
                                                          deviceField->ex2, deviceField->ey2, deviceField->ez2);

        CHECK(cudaEventRecord(eventH, streamH));

        CHECK(cudaMemcpyAsync(field->hx, deviceField->hx, bytesCount, cudaMemcpyHostToDevice, streamH));
        CHECK(cudaMemcpyAsync(field->hy, deviceField->hy, bytesCount, cudaMemcpyHostToDevice, streamH));
        CHECK(cudaMemcpyAsync(field->hz, deviceField->hz, bytesCount, cudaMemcpyHostToDevice, streamH));

        // Spawn copy thread
        hCopyParams = (CopyParams *)malloc(sizeof(CopyParams));
        hCopyParams->xSource = field->hx;
        hCopyParams->ySource = field->hy;
        hCopyParams->zSource = field->hz;
        hCopyParams->params = params;
        hCopyParams->stream = streamH;

        hThread = (pthread_t *)malloc(sizeof(pthread_t));
        pthread_create(hThread, NULL, copyResultsWithParams, hCopyParams);

        // D field
        CHECK(cudaStreamWaitEvent(streamD, eventH, 0));

        if(dThread != NULL)
            pthread_join(*dThread, NULL);

        updateDField<<<gridSize, blockSize, 0, streamD>>>(deviceField->dx0, deviceField->dy0, deviceField->dz0, 
                                                          deviceField->dx2, deviceField->dy2, deviceField->dz2, 
                                                          deviceField->hx,  deviceField->hy,  deviceField->hz);
 
        updateSources<<<gridSize, blockSize, 0, streamD>>>(deviceField->dz0, deviceField->dz2,
                                                           deviceField->hx,  deviceField->hy,
                                                           i);

        CHECK(cudaEventRecord(eventD, streamD));

        CHECK(cudaMemcpyAsync(field->dx0, deviceField->dx0, bytesCount, cudaMemcpyDeviceToHost, streamD))
        CHECK(cudaMemcpyAsync(field->dy0, deviceField->dy0, bytesCount, cudaMemcpyDeviceToHost, streamD))
        CHECK(cudaMemcpyAsync(field->dz0, deviceField->dz0, bytesCount, cudaMemcpyDeviceToHost, streamD))

        // Spawn copy thread
        dCopyParams = (CopyParams *)malloc(sizeof(CopyParams));
        dCopyParams->xSource = field->dx0;
        dCopyParams->ySource = field->dy0;
        dCopyParams->zSource = field->dz0;
        dCopyParams->params = params;
        dCopyParams->stream = streamD;

        dThread = (pthread_t *)malloc(sizeof(pthread_t));
        pthread_create(dThread, NULL, copyResultsWithParams, dCopyParams);

        // E field
        CHECK(cudaStreamWaitEvent(streamE, eventD, 0));

        if(eThread != NULL)
            pthread_join(*eThread, NULL);

        updateEField<<<gridSize, blockSize, 0, streamE>>>(deviceField->ex0, deviceField->ey0, deviceField->ez0, 
                                                          deviceField->ex2, deviceField->ey2, deviceField->ez2, 
                                                          deviceField->ex1, deviceField->ey1, deviceField->ez1, 
                                                          deviceField->dx0, deviceField->dy0, deviceField->dz0, 
                                                          deviceField->dx2, deviceField->dy2, deviceField->dz2, 
                                                          deviceField->dx1, deviceField->dy1, deviceField->dz1, 
                                                          deviceField->sigma, deviceField->epsI, deviceField->epsS, deviceField->tauD);
            
        updateMurBoundary<<<gridSize, blockSize, 0, streamE>>>(deviceField->ex0,  deviceField->ey0,  deviceField->ez0,                 
                                                               deviceField->ex2,  deviceField->ey2,  deviceField->ez2,                 
                                                               deviceField->rpx0, deviceField->rpy0, deviceField->rpz0,                         
                                                               deviceField->rpxEnd, deviceField->rpyEnd, deviceField->rpzEnd);

        CHECK(cudaEventRecord(eventE, streamE));

        CHECK(cudaMemcpyAsync(field->ex0, deviceField->ex0, bytesCount, cudaMemcpyDeviceToHost, streamE))
        CHECK(cudaMemcpyAsync(field->ey0, deviceField->ey0, bytesCount, cudaMemcpyDeviceToHost, streamE))
        CHECK(cudaMemcpyAsync(field->ez0, deviceField->ez0, bytesCount, cudaMemcpyDeviceToHost, streamE))

        // Spawn copy thread
        eCopyParams = (CopyParams *)malloc(sizeof(CopyParams));
        eCopyParams->xSource = field->ex0;
        eCopyParams->ySource = field->ey0;
        eCopyParams->zSource = field->ez0;
        eCopyParams->params = params;
        eCopyParams->stream = streamD;

        eThread = (pthread_t *)malloc(sizeof(pthread_t));
        pthread_create(eThread, NULL, copyResultsWithParams, eCopyParams);

        //Spawn write results thread
        resultsParams = (ResultsParams *)malloc(sizeof(ResultsParams));
        resultsParams->params = params;
        resultsParams->field = field;
        resultsParams->hParams = hCopyParams;
        resultsParams->dParams = dCopyParams;
        resultsParams->eParams = eCopyParams;
        resultsParams->hThread = hThread;
        resultsParams->dThread = dThread;
        resultsParams->eThread = eThread;
        resultsParams->currentIteration = i;

        pthread_create(&threads[i], NULL, writeResultsWithParams, resultsParams);

        // Run 1
        printf("Running iteration %d\n", i+1);

        // H field
        CHECK(cudaStreamWaitEvent(streamH, eventE, 0));

        if(hThread != NULL)
            pthread_join(*hThread, NULL);

        updateHField<<<gridSize, blockSize, 0, streamH>>>(deviceField->hx,  deviceField->hy,  deviceField->hz,                    
                                                          deviceField->ex0, deviceField->ey0, deviceField->ez0);

        CHECK(cudaEventRecord(eventH, streamH));
    
        CHECK(cudaMemcpyAsync(field->hx, deviceField->hx, bytesCount, cudaMemcpyHostToDevice, streamH));
        CHECK(cudaMemcpyAsync(field->hy, deviceField->hy, bytesCount, cudaMemcpyHostToDevice, streamH));
        CHECK(cudaMemcpyAsync(field->hz, deviceField->hz, bytesCount, cudaMemcpyHostToDevice, streamH));

        // Spawn copy thread
        hCopyParams = (CopyParams *)malloc(sizeof(CopyParams));
        hCopyParams->xSource = field->hx;
        hCopyParams->ySource = field->hy;
        hCopyParams->zSource = field->hz;
        hCopyParams->params = params;
        hCopyParams->stream = streamH;

        hThread = (pthread_t *)malloc(sizeof(pthread_t));
        pthread_create(hThread, NULL, copyResultsWithParams, hCopyParams);

        // D field
        CHECK(cudaStreamWaitEvent(streamD, eventH, 0));

        updateDField<<<gridSize, blockSize, 0, streamD>>>(deviceField->dx1, deviceField->dy1, deviceField->dz1, 
                                                          deviceField->dx0, deviceField->dy0, deviceField->dz0, 
                                                          deviceField->hx,  deviceField->hy,  deviceField->hz);
 
        updateSources<<<gridSize, blockSize, 0, streamD>>>(deviceField->dz1, deviceField->dz0,
                                                           deviceField->hx,  deviceField->hy,
                                                           i);

        CHECK(cudaEventRecord(eventD, streamD));

        CHECK(cudaMemcpyAsync(field->dx0, deviceField->dx1, bytesCount, cudaMemcpyDeviceToHost, streamD))
        CHECK(cudaMemcpyAsync(field->dy0, deviceField->dy1, bytesCount, cudaMemcpyDeviceToHost, streamD))
        CHECK(cudaMemcpyAsync(field->dz0, deviceField->dz1, bytesCount, cudaMemcpyDeviceToHost, streamD))

        // Spawn copy thread
        dCopyParams = (CopyParams *)malloc(sizeof(CopyParams));
        dCopyParams->xSource = field->dx0;
        dCopyParams->ySource = field->dy0;
        dCopyParams->zSource = field->dz0;
        dCopyParams->params = params;
        dCopyParams->stream = streamD;

        dThread = (pthread_t *)malloc(sizeof(pthread_t));
        pthread_create(dThread, NULL, copyResultsWithParams, dCopyParams);
 
        // E field
        CHECK(cudaStreamWaitEvent(streamE, eventD, 0));

        updateEField<<<gridSize, blockSize, 0, streamE>>>(deviceField->ex1, deviceField->ey1, deviceField->ez1,
                                                          deviceField->ex0, deviceField->ey0, deviceField->ez0,
                                                          deviceField->ex2, deviceField->ey2, deviceField->ez2,
                                                          deviceField->dx1, deviceField->dy1, deviceField->dz1,
                                                          deviceField->dx0, deviceField->dy0, deviceField->dz0,
                                                          deviceField->dx2, deviceField->dy2, deviceField->dz2,
                                                          deviceField->sigma, deviceField->epsI, deviceField->epsS, deviceField->tauD);
            
        updateMurBoundary<<<gridSize, blockSize, 0, streamE>>>(deviceField->ex1,  deviceField->ey1,  deviceField->ez1,                 
                                                               deviceField->ex0,  deviceField->ey0,  deviceField->ez0,                 
                                                               deviceField->rpx0, deviceField->rpy0, deviceField->rpz0,                         
                                                               deviceField->rpxEnd, deviceField->rpyEnd, deviceField->rpzEnd);

        CHECK(cudaEventRecord(eventE, streamE));

        CHECK(cudaMemcpyAsync(field->ex0, deviceField->ex1, bytesCount, cudaMemcpyDeviceToHost, streamE))
        CHECK(cudaMemcpyAsync(field->ey0, deviceField->ey1, bytesCount, cudaMemcpyDeviceToHost, streamE))
        CHECK(cudaMemcpyAsync(field->ez0, deviceField->ez1, bytesCount, cudaMemcpyDeviceToHost, streamE))

        // Spawn copy thread
        eCopyParams = (CopyParams *)malloc(sizeof(CopyParams));
        eCopyParams->xSource = field->ex0;
        eCopyParams->ySource = field->ey0;
        eCopyParams->zSource = field->ez0;
        eCopyParams->params = params;
        eCopyParams->stream = streamD;

        eThread = (pthread_t *)malloc(sizeof(pthread_t));
        pthread_create(eThread, NULL, copyResultsWithParams, eCopyParams);

        //Spawn write results thread
        resultsParams = (ResultsParams *)malloc(sizeof(ResultsParams));
        resultsParams->params = params;
        resultsParams->field = field;
        resultsParams->hParams = hCopyParams;
        resultsParams->dParams = dCopyParams;
        resultsParams->eParams = eCopyParams;
        resultsParams->hThread = hThread;
        resultsParams->dThread = dThread;
        resultsParams->eThread = eThread;
        resultsParams->currentIteration = i+1;

        pthread_create(&threads[i+1], NULL, writeResultsWithParams, resultsParams);

        // Run 2
        printf("Running iteration %d\n", i+2);

        // H field
        CHECK(cudaStreamWaitEvent(streamH, eventE, 0));

        if(hThread != NULL)
            pthread_join(*hThread, NULL);

        updateHField<<<gridSize, blockSize, 0, streamH>>>(deviceField->hx,  deviceField->hy,  deviceField->hz,                    
                                                          deviceField->ex1, deviceField->ey1, deviceField->ez1);

        CHECK(cudaEventRecord(eventH, streamH));

        CHECK(cudaMemcpyAsync(field->hx, deviceField->hx, bytesCount, cudaMemcpyHostToDevice, streamH));
        CHECK(cudaMemcpyAsync(field->hy, deviceField->hy, bytesCount, cudaMemcpyHostToDevice, streamH));
        CHECK(cudaMemcpyAsync(field->hz, deviceField->hz, bytesCount, cudaMemcpyHostToDevice, streamH));

        // Spawn copy thread
        hCopyParams = (CopyParams *)malloc(sizeof(CopyParams));
        hCopyParams->xSource = field->hx;
        hCopyParams->ySource = field->hy;
        hCopyParams->zSource = field->hz;
        hCopyParams->params = params;
        hCopyParams->stream = streamH;

        hThread = (pthread_t *)malloc(sizeof(pthread_t));
        pthread_create(hThread, NULL, copyResultsWithParams, hCopyParams);

        // D field
        CHECK(cudaStreamWaitEvent(streamD, eventH, 0));

        if(dThread != NULL)
            pthread_join(*dThread, NULL);

        updateDField<<<gridSize, blockSize, 0, streamD>>>(deviceField->dx2, deviceField->dy2, deviceField->dz2, 
                                                          deviceField->dx1, deviceField->dy1, deviceField->dz1, 
                                                          deviceField->hx,  deviceField->hy,  deviceField->hz);
 
        updateSources<<<gridSize, blockSize, 0, streamD>>>(deviceField->dz2, deviceField->dz1,
                                                           deviceField->hx,  deviceField->hy,
                                                           i);

        CHECK(cudaEventRecord(eventD, streamD));

        CHECK(cudaMemcpyAsync(field->dx0, deviceField->dx2, bytesCount, cudaMemcpyDeviceToHost, streamD))
        CHECK(cudaMemcpyAsync(field->dy0, deviceField->dy2, bytesCount, cudaMemcpyDeviceToHost, streamD))
        CHECK(cudaMemcpyAsync(field->dz0, deviceField->dz2, bytesCount, cudaMemcpyDeviceToHost, streamD))

        // Spawn copy thread
        dCopyParams = (CopyParams *)malloc(sizeof(CopyParams));
        dCopyParams->xSource = field->dx0;
        dCopyParams->ySource = field->dy0;
        dCopyParams->zSource = field->dz0;
        dCopyParams->params = params;
        dCopyParams->stream = streamD;

        dThread = (pthread_t *)malloc(sizeof(pthread_t));
        pthread_create(dThread, NULL, copyResultsWithParams, dCopyParams);
            
        // E field
        CHECK(cudaStreamWaitEvent(streamE, eventD, 0));

        if(eThread != NULL)
            pthread_join(*eThread, NULL);

        updateEField<<<gridSize, blockSize, 0, streamE>>>(deviceField->ex2, deviceField->ey2, deviceField->ez2, 
                                                          deviceField->ex1, deviceField->ey1, deviceField->ez1, 
                                                          deviceField->ex0, deviceField->ey0, deviceField->ez0, 
                                                          deviceField->dx2, deviceField->dy2, deviceField->dz2, 
                                                          deviceField->dx1, deviceField->dy1, deviceField->dz1, 
                                                          deviceField->dx0, deviceField->dy0, deviceField->dz0, 
                                                          deviceField->sigma, deviceField->epsI, deviceField->epsS, deviceField->tauD);
            
        updateMurBoundary<<<gridSize, blockSize, 0, streamE>>>(deviceField->ex2,  deviceField->ey2,  deviceField->ez2,                 
                                                               deviceField->ex1,  deviceField->ey1,  deviceField->ez1,                 
                                                               deviceField->rpx0, deviceField->rpy0, deviceField->rpz0,                         
                                                               deviceField->rpxEnd, deviceField->rpyEnd, deviceField->rpzEnd);

        CHECK(cudaEventRecord(eventE, streamE));

        CHECK(cudaMemcpyAsync(field->ex0, deviceField->ex2, bytesCount, cudaMemcpyDeviceToHost, streamE))
        CHECK(cudaMemcpyAsync(field->ey0, deviceField->ey2, bytesCount, cudaMemcpyDeviceToHost, streamE))
        CHECK(cudaMemcpyAsync(field->ez0, deviceField->ez2, bytesCount, cudaMemcpyDeviceToHost, streamE))

        // Spawn copy thread
        eCopyParams = (CopyParams *)malloc(sizeof(CopyParams));
        eCopyParams->xSource = field->ex0;
        eCopyParams->ySource = field->ey0;
        eCopyParams->zSource = field->ez0;
        eCopyParams->params = params;
        eCopyParams->stream = streamD;

        eThread = (pthread_t *)malloc(sizeof(pthread_t));
        pthread_create(eThread, NULL, copyResultsWithParams, eCopyParams);

        //Spawn write results thread
        resultsParams = (ResultsParams *)malloc(sizeof(ResultsParams));
        resultsParams->params = params;
        resultsParams->field = field;
        resultsParams->hParams = hCopyParams;
        resultsParams->dParams = dCopyParams;
        resultsParams->eParams = eCopyParams;
        resultsParams->hThread = hThread;
        resultsParams->dThread = dThread;
        resultsParams->eThread = eThread;
        resultsParams->currentIteration = i+2;

        pthread_create(&threads[i+2], NULL, writeResultsWithParams, resultsParams);
    }


    // Wait for all threads to finish
    for(int i=0; i<params->iterationsCount; i++) {
        pthread_join(threads[i], NULL);
    }

    // Clean up
    //free(threads);

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
                                            &params->sources[i*3 + 0],
                                            &params->sources[i*3 + 1],
                                            &params->sources[i*3 + 2]);
        params->sources[i*3 + 0] -= 1;
        params->sources[i*3 + 1] -= 1;
        params->sources[i*3 + 2] -= 1;
    }
    //unused (pulse_type)
    fgets(temp, tempLength, paramsFile);
    //fsigma (sigma)
    fscanf(paramsFile, "%s %f\n", temp, &params->defaultSigma);
    //feps_s (eps_s)
    fscanf(paramsFile, "%s %f\n", temp, &params->defaultEpsS);
    //feps_inf (eps_i)
    fscanf(paramsFile, "%s %f\n", temp, &params->defaultEpsI);
    //ftau_d (tau_d)
    fscanf(paramsFile, "%s %f\n", temp, &params->defaultTauD);
    
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
    params->eps0 = 1.0 / (params->mu0 * params->c * params->c);

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
    printf("Field size:                 %04dx%04dx%04d\n", params->nx, params->ny, params->nz);
    printf("Iterations count:           %d\n", params->iterationsCount);
    printf("Input path:                 %s\n", params->inputPath);
    printf("Output path:                %s\n", params->outputPath);
    printf("Elements per wavelength:    %d\n", params->elementsPerWave);
    printf("Wave frequency:             %9.3E\n", params->waveFrequency);
    printf("Pulse width:                %9.3E\n", params->pulseWidth);
    printf("Pulse modulation frequency: %9.3E\n", params->pulseModulationFrequency);
    printf("Sources count:              %d\n", params->sourcesCount);
    for(int i=0; i<params->sourcesCount; i++)
        printf("Source position:            %04dx%04dx%04d\n", params->sources[i*3 + 0] + 1,
                                                               params->sources[i*3 + 1] + 1,
                                                               params->sources[i*3 + 2] + 1);
    printf("Default sigma:              %9.3E\n", params->defaultSigma);
    printf("Default eps_s:              %9.3E\n", params->defaultEpsS);
    printf("Default eps_i:              %9.3E\n", params->defaultEpsI);
    printf("Default tau_d:              %9.3E\n", params->defaultTauD);
    printf("\n");

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    printf("Running on %s\n", deviceProp.name);
    printf("Compute capability: %d.%d", deviceProp.major, deviceProp.minor);
    printf("Memory available: %.2f", deviceProp.totalGlobalMem / (1024.0 * 1024.0));

    int usedBytes = params->nx * params->ny * params->nz;
    usedBytes *= 7 * 3 + 10; // e1, e2, e3, h, d0, d1, d2, d3, eps, tau, sigma, rp, etc...
    usedBytes *= sizeof(float);
    printf("Memory requirements: %.2f MB", (float)usedBytes / (1024.0 * 1024.0));
}


FdtdField *initFieldWithParams(FdtdParams *params)
{
    int n = params->nx * params->ny * params->nz; 

    FdtdField *field = (FdtdField *)malloc(sizeof(FdtdField));
    if(field == NULL) {
        printf("Couldn't allocate field\n");
        exit(EXIT_FAILURE);
    }

    //H
    CHECK(cudaHostAlloc(&field->hx, n * sizeof(float), cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->hy, n * sizeof(float), cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->hz, n * sizeof(float), cudaHostAllocDefault))

    //D
    CHECK(cudaHostAlloc(&field->dx0, n * sizeof(float), cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dy0, n * sizeof(float), cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->dz0, n * sizeof(float), cudaHostAllocDefault))

    //E
    CHECK(cudaHostAlloc(&field->ex0, n * sizeof(float), cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ey0, n * sizeof(float), cudaHostAllocDefault))
    CHECK(cudaHostAlloc(&field->ez0, n * sizeof(float), cudaHostAllocDefault))

    // sigma, eps, tau
    //CHECK(cudaHostAlloc(&field->sigma, n * sizeof(float), cudaHostAllocDefault))
    //CHECK(cudaHostAlloc(&field->epsS,  n * sizeof(float), cudaHostAllocDefault))
    //CHECK(cudaHostAlloc(&field->epsI,  n * sizeof(float), cudaHostAllocDefault))
    //CHECK(cudaHostAlloc(&field->tauD,  n * sizeof(float), cudaHostAllocDefault))
    
    field->sigma = (float *)malloc( n * sizeof(float));
    field->epsS  = (float *)malloc( n * sizeof(float));
    field->epsI  = (float *)malloc( n * sizeof(float));
    field->tauD  = (float *)malloc( n * sizeof(float));

    for(int i = 0; i < n; i++) {
        field->sigma[i] = params->defaultSigma;
        field->epsS[i]  = params->defaultEpsS;
        field->epsI[i]  = params->defaultEpsI;
        field->tauD[i]  = params->defaultTauD;
    }

    // rp
    //CHECK(cudaHostAlloc(&field->rpx0, n * sizeof(float), cudaHostAllocDefault))
    //CHECK(cudaHostAlloc(&field->rpy0, n * sizeof(float), cudaHostAllocDefault))
    //CHECK(cudaHostAlloc(&field->rpz0, n * sizeof(float), cudaHostAllocDefault))

    //CHECK(cudaHostAlloc(&field->rpxEnd, n * sizeof(float), cudaHostAllocDefault))
    //CHECK(cudaHostAlloc(&field->rpyEnd, n * sizeof(float), cudaHostAllocDefault))
    //CHECK(cudaHostAlloc(&field->rpzEnd, n * sizeof(float), cudaHostAllocDefault))

    field->rpx0 = (float *)malloc(n * sizeof(float));
    field->rpy0 = (float *)malloc(n * sizeof(float)); 
    field->rpz0 = (float *)malloc(n * sizeof(float)); 

    field->rpxEnd = (float *)malloc(n * sizeof(float));
    field->rpyEnd = (float *)malloc(n * sizeof(float));
    field->rpzEnd = (float *)malloc(n * sizeof(float));  

    return field;
}


void deallocField(FdtdField *field)
{
    //H
    CHECK(cudaFreeHost(field->hx));
    CHECK(cudaFreeHost(field->hy));
    CHECK(cudaFreeHost(field->hz));

    //D
    CHECK(cudaFreeHost(field->dx0));
    CHECK(cudaFreeHost(field->dy0));
    CHECK(cudaFreeHost(field->dz0));

    //E
    CHECK(cudaFreeHost(field->ex0));
    CHECK(cudaFreeHost(field->ey0));
    CHECK(cudaFreeHost(field->ez0));

    //sigma, eps, tau
    free(field->sigma);
    free(field->epsS);
    free(field->epsI);
    free(field->tauD);

    //rp
    free(field->rpx0);
    free(field->rpy0);
    free(field->rpz0);

    free(field->rpxEnd);
    free(field->rpyEnd);
    free(field->rpzEnd);

    free(field);
}


FdtdField *initDeviceFieldWithParams(FdtdParams *params)
{
    int n = params->nx * params->ny * params->nz; 

    FdtdField *field = (FdtdField *)malloc(sizeof(FdtdField));

    // e
    CHECK(cudaMalloc(&field->ex0, n * sizeof(float)))
    CHECK(cudaMalloc(&field->ey0, n * sizeof(float)))
    CHECK(cudaMalloc(&field->ez0, n * sizeof(float)))

    CHECK(cudaMalloc(&field->ex1, n * sizeof(float)))
    CHECK(cudaMalloc(&field->ey1, n * sizeof(float)))
    CHECK(cudaMalloc(&field->ez1, n * sizeof(float)))

    CHECK(cudaMalloc(&field->ex2, n * sizeof(float)))
    CHECK(cudaMalloc(&field->ey2, n * sizeof(float)))
    CHECK(cudaMalloc(&field->ez2, n * sizeof(float)))

    // h
    CHECK(cudaMalloc(&field->hx, n * sizeof(float)))
    CHECK(cudaMalloc(&field->hy, n * sizeof(float)))
    CHECK(cudaMalloc(&field->hz, n * sizeof(float)))

    // d
    CHECK(cudaMalloc(&field->dx0, n * sizeof(float)))
    CHECK(cudaMalloc(&field->dy0, n * sizeof(float)))
    CHECK(cudaMalloc(&field->dz0, n * sizeof(float)))

    CHECK(cudaMalloc(&field->dx1, n * sizeof(float)))
    CHECK(cudaMalloc(&field->dy1, n * sizeof(float)))
    CHECK(cudaMalloc(&field->dz1, n * sizeof(float)))

    CHECK(cudaMalloc(&field->dx2, n * sizeof(float)))
    CHECK(cudaMalloc(&field->dy2, n * sizeof(float)))
    CHECK(cudaMalloc(&field->dz2, n * sizeof(float)))

    // sigma, eps, tau
    CHECK(cudaMalloc(&field->epsI,  n * sizeof(float)))
    CHECK(cudaMalloc(&field->epsS,  n * sizeof(float)))
    CHECK(cudaMalloc(&field->tauD,  n * sizeof(float)))
    CHECK(cudaMalloc(&field->sigma, n * sizeof(float)))

    // rp
    CHECK(cudaMalloc(&field->rpx0, n * sizeof(float)))
    CHECK(cudaMalloc(&field->rpy0, n * sizeof(float)))
    CHECK(cudaMalloc(&field->rpz0, n * sizeof(float)))

    CHECK(cudaMalloc(&field->rpxEnd, n * sizeof(float)))
    CHECK(cudaMalloc(&field->rpyEnd, n * sizeof(float)))
    CHECK(cudaMalloc(&field->rpzEnd, n * sizeof(float)))

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
        printf("Couldn't alocate %ld bytes for specs\n", (long)specsCount*4*sizeof(float));
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

    for(int i=0; i<specsCount; i++) {
        fscanf(specsFile, "%d %s %g %g %g %g\n", &index, temp, &sigmaValue, &epsSValue, &epsIValue, &tauDValue);

        specs[index*4 + 0] = sigmaValue;
        specs[index*4 + 1] = epsSValue;
        specs[index*4 + 2] = epsIValue;
        specs[index*4 + 3] = tauDValue;

        if(index >= specsCount)
            break;
    }

    //fclose(specsFile);

    // Load materials
    for(int iz=0; iz<params->nz; iz++) {
        char materialFileName[1024];
        sprintf(materialFileName, "%s/v1_%05d.pgm", materialsPath, iz+1);
        FILE *materialFile = fopen(materialFileName, "r");
        
        if(materialFile == NULL) {
            printf("Couldn\'t open file %s\n", materialFileName);
            exit(EXIT_FAILURE);
        }

        int width, height;
        fscanf(materialFile, "%s %s %s %d %d %s", temp, temp, temp, &width, &height, temp);

        for(int iy=0; iy<params->ny; iy++) {
            for(int ix=0; ix<params->nx; ix++) {
                int code;
                fscanf(materialFile, "%d", &code);

                int offset = iz*params->nx*params->ny + iy*params->nx + ix;
                field->sigma[offset] = specs[code*4 + 0];
                field->epsS[offset] =  specs[code*4 + 1];
                field->epsI[offset] =  specs[code*4 + 2];
                field->tauD[offset] =  specs[code*4 + 3];
            }
        }

        fclose(materialFile);
    }

    //free(specs);
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


void setupSources(FdtdParams *params)
{
    int fine, temp, i2, istart;
    float *tmpdata, *tmpdata2;
    int tmpOff = 1<<16;

    params->jzCount = tmpOff;
    params->jz = (float *)calloc(tmpOff,     sizeof(float));
    tmpdata    = (float *)calloc(tmpOff * 2, sizeof(float));
    tmpdata2   = (float *)calloc(tmpOff * 2, sizeof(float));
    
    //fine & temp
    fine = (1<<13) * params->pulseWidth * params->waveFrequency * params->dt;
    temp = 1.0/(params->pulseWidth * params->waveFrequency)/(params->dt / fine)/2.0;
    
    //tmpdata
    for(int i = -temp - 1; i <= temp + 1; i++) {
        float v1 = ((float)i/(((float)temp + 1.0)/4.0));
        float v2 = exp(-pow(v1, 2.0));
        float v3 = cos(2.0 * acos(-1.0) * params->pulseModulationFrequency * params->waveFrequency * i * (params->dt / fine));
        tmpdata[i + tmpOff] = v2 * v3;
    }

    //istart
    for(int i = -(1<<12); i < (1<<12); i++) {
         if((fabs(tmpdata[i + tmpOff]) >= 1e-9) && (i % fine == 0)) {
            istart = i;
            break;
         }
    }
    
    //setup jz 1/2
    i2 = 0;
    for(int i = istart; i <= temp+1; i += fine) {
        float val = tmpdata[i + tmpOff] * 1e-15 / params->dt / 3.0;
        params->jz[i2] = val;
        i2++;
    }
    
    //setup tmpdata2
    for(int i = 2; i <= (1<<14); i++) {
        float val = (((params->jz[i + 1 - 1] - params->jz[i - 1]) / params->dt) +
                     ((params->jz[i - 1] - params->jz[i - 1 - 1]) / params->dt)) / 
                    2.0 * (params->dt * params->dz) / (params->dx * params->dy * params->dz);
                                    
        tmpdata2[i - 1 + tmpOff] = val;
    }
    
    //setup jz 2/2
    for(int i=0; i < 1<<14; i++) {
        params->jz[i] = tmpdata2[i + tmpOff + 1];
    }

    free(tmpdata2);
    free(tmpdata);
}


void copyDataToDevice(FdtdParams *params, FdtdField *field, FdtdField *deviceField)
{
    int bytesCount = params->nx * params->ny * params->nz * sizeof(float); 

    //H
    CHECK(cudaMemset(deviceField->hx, 0, bytesCount))
    CHECK(cudaMemset(deviceField->hy, 0, bytesCount))
    CHECK(cudaMemset(deviceField->hz, 0, bytesCount))

    //D
    CHECK(cudaMemset(deviceField->dx0, 0, bytesCount))
    CHECK(cudaMemset(deviceField->dy0, 0, bytesCount))
    CHECK(cudaMemset(deviceField->dz0, 0, bytesCount))

    CHECK(cudaMemset(deviceField->dx1, 0, bytesCount))
    CHECK(cudaMemset(deviceField->dy1, 0, bytesCount))
    CHECK(cudaMemset(deviceField->dz1, 0, bytesCount))

    CHECK(cudaMemset(deviceField->dx2, 0, bytesCount))
    CHECK(cudaMemset(deviceField->dy2, 0, bytesCount))
    CHECK(cudaMemset(deviceField->dz2, 0, bytesCount))

    //E
    CHECK(cudaMemset(deviceField->ex0, 0, bytesCount))
    CHECK(cudaMemset(deviceField->ey0, 0, bytesCount))
    CHECK(cudaMemset(deviceField->ez0, 0, bytesCount))

    CHECK(cudaMemset(deviceField->ex1, 0, bytesCount))
    CHECK(cudaMemset(deviceField->ey1, 0, bytesCount))
    CHECK(cudaMemset(deviceField->ez1, 0, bytesCount))

    CHECK(cudaMemset(deviceField->ex2, 0, bytesCount))
    CHECK(cudaMemset(deviceField->ey2, 0, bytesCount))
    CHECK(cudaMemset(deviceField->ez2, 0, bytesCount))

    //eps, tau, sigma
    CHECK(cudaMemcpy(deviceField->epsI,  field->epsI,  bytesCount, cudaMemcpyHostToDevice))
    CHECK(cudaMemcpy(deviceField->epsS,  field->epsS,  bytesCount, cudaMemcpyHostToDevice))
    CHECK(cudaMemcpy(deviceField->tauD,  field->tauD,  bytesCount, cudaMemcpyHostToDevice))
    CHECK(cudaMemcpy(deviceField->sigma, field->sigma, bytesCount, cudaMemcpyHostToDevice))

    CHECK(cudaMemcpy(deviceField->rpx0, field->rpx0, bytesCount, cudaMemcpyHostToDevice))
    CHECK(cudaMemcpy(deviceField->rpy0, field->rpy0, bytesCount, cudaMemcpyHostToDevice))
    CHECK(cudaMemcpy(deviceField->rpz0, field->rpz0, bytesCount, cudaMemcpyHostToDevice))

    CHECK(cudaMemcpy(deviceField->rpxEnd, field->rpxEnd, bytesCount, cudaMemcpyHostToDevice))
    CHECK(cudaMemcpy(deviceField->rpyEnd, field->rpyEnd, bytesCount, cudaMemcpyHostToDevice))
    CHECK(cudaMemcpy(deviceField->rpzEnd, field->rpzEnd, bytesCount, cudaMemcpyHostToDevice))
}


void *copyResultsWithParams(void *params)
{
    CopyParams *copyParams = (CopyParams *)params;

    int bytesCount = copyParams->params->nx * copyParams->params->ny * copyParams->params->nz * sizeof(float);

    copyParams->xBuffer = (float *)malloc(bytesCount);
    copyParams->yBuffer = (float *)malloc(bytesCount);
    copyParams->zBuffer = (float *)malloc(bytesCount);

    CHECK(cudaStreamSynchronize(copyParams->stream))

    memcpy(copyParams->xBuffer, copyParams->xSource, bytesCount);
    memcpy(copyParams->yBuffer, copyParams->ySource, bytesCount);
    memcpy(copyParams->zBuffer, copyParams->zSource, bytesCount);

    pthread_exit(NULL);
}


void *writeResultsWithParams(void *params)
{
    ResultsParams *resultsParams = (ResultsParams *)params;

    pthread_join(*resultsParams->hThread, NULL);
    pthread_join(*resultsParams->dThread, NULL);
    pthread_join(*resultsParams->eThread, NULL);

    writeResults(resultsParams->params, resultsParams->field,
                 resultsParams->hParams->xBuffer, resultsParams->hParams->yBuffer, resultsParams->hParams->zBuffer,
                 resultsParams->dParams->xBuffer, resultsParams->dParams->yBuffer, resultsParams->dParams->zBuffer,
                 resultsParams->eParams->xBuffer, resultsParams->eParams->yBuffer, resultsParams->eParams->zBuffer,
                 resultsParams->currentIteration);

    free(resultsParams->hParams);
    free(resultsParams->dParams);
    free(resultsParams->eParams);

    free(params);

    pthread_exit(NULL);
}


void writeResults(FdtdParams *params, FdtdField *field,
                  float *hxSource, float *hySource, float *hzSource,
                  float *dxSource, float *dySource, float *dzSource,
                  float *exSource, float *eySource, float *ezSource,
                  int currentIteration)
{
    char outputFilePath[1024];
    FILE *outputFile;

    // Used by OFFSET macro
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    // Output hx
    for(int iz = 0; iz < params->nz; iz++) {
        sprintf(outputFilePath, "%s/H_field_x_%05d_z%05d.out", params->outputPath, currentIteration, iz);

        outputFile = fopen(outputFilePath, "w");
        if(outputFile == NULL) {
            printf("Couldn\'t open file %s\n", outputFilePath);
            exit(EXIT_FAILURE);
        }

        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix=0; ix < params->nx; ix++)
                fprintf(outputFile, "% 9.3E", OFFSET(hxSource, ix, iy, iz));
            
            fprintf(outputFile, "\n");
        }

        fclose(outputFile);
    }

    // Output hy
    for(int iz = 0; iz < params->nz; iz++) {
        sprintf(outputFilePath, "%s/H_field_y_%05d_z%05d.out", params->outputPath, currentIteration, iz);

        outputFile = fopen(outputFilePath, "w");
        if(outputFile == NULL) {
            printf("Couldn\'t open file %s\n", outputFilePath);
            exit(EXIT_FAILURE);
        }

        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix=0; ix < params->nx; ix++)
                fprintf(outputFile, "% 9.3E", OFFSET(hySource, ix, iy, iz));
            
            fprintf(outputFile, "\n");
        }

        fclose(outputFile);
    }

    // Output hz
    for(int iz = 0; iz < params->nz; iz++) {
        sprintf(outputFilePath, "%s/H_field_z_%05d_z%05d.out", params->outputPath, currentIteration, iz);

        outputFile = fopen(outputFilePath, "w");
        if(outputFile == NULL) {
            printf("Couldn\'t open file %s\n", outputFilePath);
            exit(EXIT_FAILURE);
        }

        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix=0; ix < params->nx; ix++)
                fprintf(outputFile, "% 9.3E", OFFSET(hzSource, ix, iy, iz));
            
            fprintf(outputFile, "\n");
        }

        fclose(outputFile);
    }

    // Output dx
    for(int iz = 0; iz < params->nz; iz++) {
        sprintf(outputFilePath, "%s/D_field_x_%05d_z%05d.out", params->outputPath, currentIteration, iz);

        outputFile = fopen(outputFilePath, "w");
        if(outputFile == NULL) {
            printf("Couldn\'t open file %s\n", outputFilePath);
            exit(EXIT_FAILURE);
        }

        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix=0; ix < params->nx; ix++)
                fprintf(outputFile, "% 9.3E", OFFSET(dxSource, ix, iy, iz));
            
            fprintf(outputFile, "\n");
        }

        fclose(outputFile);
    }

    // Output dy
    for(int iz = 0; iz < params->nz; iz++) {
        sprintf(outputFilePath, "%s/D_field_y_%05d_z%05d.out", params->outputPath, currentIteration, iz);

        outputFile = fopen(outputFilePath, "w");
        if(outputFile == NULL) {
            printf("Couldn\'t open file %s\n", outputFilePath);
            exit(EXIT_FAILURE);
        }

        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix=0; ix < params->nx; ix++)
                fprintf(outputFile, "% 9.3E", OFFSET(dySource, ix, iy, iz));
            
            fprintf(outputFile, "\n");
        }

        fclose(outputFile);
    }

    // Output dz
    for(int iz = 0; iz < params->nz; iz++) {
        sprintf(outputFilePath, "%s/D_field_z_%05d_z%05d.out", params->outputPath, currentIteration, iz);

        outputFile = fopen(outputFilePath, "w");
        if(outputFile == NULL) {
            printf("Couldn\'t open file %s\n", outputFilePath);
            exit(EXIT_FAILURE);
        }

        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix=0; ix < params->nx; ix++)
                fprintf(outputFile, "% 9.3E", OFFSET(dzSource, ix, iy, iz));
            
            fprintf(outputFile, "\n");
        }

        fclose(outputFile);
    }

    // Output ex
    for(int iz = 0; iz < params->nz; iz++) {
        sprintf(outputFilePath, "%s/E_field_x_%05d_z%05d.out", params->outputPath, currentIteration, iz);

        outputFile = fopen(outputFilePath, "w");
        if(outputFile == NULL) {
            printf("Couldn\'t open file %s\n", outputFilePath);
            exit(EXIT_FAILURE);
        }

        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix=0; ix < params->nx; ix++)
                fprintf(outputFile, "% 9.3E", OFFSET(exSource, ix, iy, iz));
            
            fprintf(outputFile, "\n");
        }

        fclose(outputFile);
    }

    // Output ey
    for(int iz = 0; iz < params->nz; iz++) {
        sprintf(outputFilePath, "%s/E_field_y_%05d_z%05d.out", params->outputPath, currentIteration, iz);

        outputFile = fopen(outputFilePath, "w");
        if(outputFile == NULL) {
            printf("Couldn\'t open file %s\n", outputFilePath);
            exit(EXIT_FAILURE);
        }

        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix=0; ix < params->nx; ix++)
                fprintf(outputFile, "% 9.3E", OFFSET(eySource, ix, iy, iz));
            
            fprintf(outputFile, "\n");
        }

        fclose(outputFile);
    }

    // Output ez
    for(int iz = 0; iz < params->nz; iz++) {
        sprintf(outputFilePath, "%s/E_field_z_%05d_z%05d.out", params->outputPath, currentIteration, iz);

        outputFile = fopen(outputFilePath, "w");
        if(outputFile == NULL) {
            printf("Couldn\'t open file %s\n", outputFilePath);
            exit(EXIT_FAILURE);
        }

        for(int iy = 0; iy < params->ny; iy++) {
            for(int ix=0; ix < params->nx; ix++)
                fprintf(outputFile, "% 9.3E", OFFSET(ezSource, ix, iy, iz));
            
            fprintf(outputFile, "\n");
        }

        fclose(outputFile);
    }


    //Cleanup unnecessary buffers
    free(hxSource);
    free(hySource);
    free(hzSource);

    free(dxSource);
    free(dySource);
    free(dzSource);

    free(exSource);
    free(eySource);
    free(ezSource);
}
