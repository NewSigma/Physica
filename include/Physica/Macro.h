/*
 * Copyright 2019-2024 WeiBo He.
 *
 * This file is part of Physica.
 *
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
#pragma once

#include <climits>
#include "Physica/Config.h"

#if PhysicaWordSize == 64
    #define PHYSICA_64BIT
#elif PhysicaWordSize == 32;
    #define PHYSICA_32BIT
#endif

/*!
 * Improve: Platform dependent, may be wrong on some spatial platforms, add a test before compile.
 * Note: Use of double_extract may cause several warnings in valgrind.
 */
#if PhysicaEndian == PhysicaBigEndian
    union float_extract {
        float value;
        struct {
            unsigned int sign : 1;
            unsigned int exp : 8;
            unsigned int high : 7;
            unsigned int low : 16;
        };
    };

    union double_extract {
        double value;
        struct {
            unsigned int sign : 1;
            unsigned int exp : 11;
            unsigned int high : 20;
            unsigned int low : 32;
        };
    };
#elif PhysicaEndian == PhysicaLittleEndian
    union float_extract {
        float value;
        struct {
            unsigned int low : 16;
            unsigned int high : 7;
            unsigned int exp : 8;
            unsigned int sign : 1;
        };
    };

    union double_extract {
        double value;
        struct {
            unsigned int low : 32;
            unsigned int high : 20;
            unsigned int exp : 11;
            unsigned int sign : 1;
        };
    };
#endif

#define PHYSICA_API __attribute__ ((visibility ("default")))

#ifndef PHYSICA_CUDA
    #define __host__
    #define __device__
#else
    #include <cuda_runtime_api.h>
#endif

namespace Physica {
    constexpr static bool IsHDF5Enabled() {
    #ifdef PHYSICA_HDF5
        return true;
    #else
        return false;
    #endif
    }

    constexpr static bool IsMKLEnabled() {
    #ifdef PHYSICA_MKL
        return true;
    #else
        return false;
    #endif
    }

    constexpr static bool IsCUDAEnabled() {
    #ifdef PHYSICA_CUDA
        return true;
    #else
        return false;
    #endif
    }
}
