#ifndef PHYSICA_DEBUGUTIL_H
#define PHYSICA_DEBUGUTIL_H

#include <cstring>
#include <cuda_runtime_api.h>

//e.g turn /home/user/Physica/Physica.cpp into Physica
#ifdef WIN32
#define FILENAME(x) (strrchr(x, '\\') ? strrchr(x, '\\') + 1 : x)
#else
#define FILENAME(x) (strrchr(x, '/') ? strrchr(x, '/') + 1 : x)
#endif
//IDEA: A better logger for both C++ and CUDA is required.
#define cuDebug(x) do { printf("[] [Debug] [%s:%d]: %s\n", __FILENAME__, __LINE__, x); } while(false)
#define cuWarning(x) do { printf("[] [Warning] [%s:%d]: %s\n", __FILENAME__, __LINE__, x); } while(false)
#define cuCritical(x) do { printf("[] [Critical] [%s:%d]: %s\n", __FILENAME__, __LINE__, x); } while(false)
#define cuFatal(x)                                                    \
    do {                                                              \
        printf("[] [Fatal] [%s:%d]: %s\n", __FILENAME__, __LINE__, x);\
        cudaDeviceReset();                                            \
        exit(EXIT_FAILURE);                                           \
    } while(false)
#define cuInfo(x) do { printf("[] [Info] [%s:%d]: %s\n", __FILENAME__, __LINE__, x); } while(false)
#define checkCudaError(x)                                                                                                       \
    do {                                                                                                                        \
        if(x) {                                                                                                                 \
            printf("[] [Fatal] [%s:%d]: CUDA error encountered! Error code: %s\n", FILENAME(__FILE__), __LINE__, cudaGetErrorName(x));\
            cudaDeviceReset();                                                                                                  \
            exit(EXIT_FAILURE);                                                                                                 \
        }                                                                                                                       \
    } while(false)

#endif