#ifndef _UTILS_H
#define _UTILS_H


#define CHECK(call)                                                \
{                                                                  \
    const cudaError_t error = call;                                \
    if(error != cudaSuccess) {                                     \
        printf("Error: %s: %d\n", __FILE__, __LINE__);               \
        printf("Code %d: %s\n", error, cudaGetErrorString(error)); \
        exit(EXIT_FAILURE);                                        \
    }                                                              \
}

#define OFFSET(p, x, y, z) (p[(z)*ny*nx + (y)*nx + (x)])
#define OFFSETRP(p, x, y, z) (p[(z)*rpny*nx + (y)*rpnx + (x)])

void check(int condition, char *message);

#endif
