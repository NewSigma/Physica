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

#include <thread>
#include "Physica/Core/Parallel/StreamPool.cuh"
#include "SequentialExecutor.h"

namespace Physica::Core {
    class CudaExecutor;

    namespace Internal {
        template<class T> class Traits;

        template<>
        class Traits<CudaExecutor> {
        public:
            constexpr static bool isCudaEnabled = true;
        };
    }
    /**
     * Single thread with cuda support
     */
    class CudaExecutor : public SequentialExecutor {
    public:
        static void wait() {
            while (StreamPool::getStream().query() != cudaSuccess)
                std::this_thread::yield();
            cudaCheck(cudaGetLastError());
        }
    };
}
