/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef PHYSICA_DEBUGUTIL_H
#define PHYSICA_DEBUGUTIL_H

#include <cstring>
#include <cuda_runtime_api.h>

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
            printf("[] [Fatal] [%s:%d]: CUDA error encountered! Error code: %s\n", __FILENAME__, __LINE__, cudaGetErrorName(x));\
            cudaDeviceReset();                                                                                                  \
            exit(EXIT_FAILURE);                                                                                                 \
        }                                                                                                                       \
    } while(false)

#endif