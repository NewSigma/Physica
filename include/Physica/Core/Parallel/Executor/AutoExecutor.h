/*
 * Copyright 2023 WeiBo He.
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

#include "ThreadExecutor.h"
#ifdef PHYSICA_CUDA
    #include "CudaExecutor.cuh"
#endif

namespace Physica::Core {
    class AutoExecutor;

    namespace Internal {
        template<class T> class Traits;

        template<>
        class Traits<AutoExecutor> {
        public:
            constexpr static bool isCPUEnabled = true;
        #ifdef PHYSICA_CUDA
            constexpr static bool isCudaEnabled = true;
        #else
            constexpr static bool isCudaEnabled = false;
        #endif
        };
    }
    /**
     * \class AutoExecutor is devoted to provide a load balancing heterogeneous computing.
     */
    class AutoExecutor : public ThreadExecutor {
        using Base = ThreadExecutor;
    public:
        static void wait() {
        #ifdef PHYSICA_CUDA
            CudaExecutor::wait();
        #endif
        }
    };
}
