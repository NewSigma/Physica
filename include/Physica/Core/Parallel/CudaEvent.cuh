/*
 * Copyright 2023-2024 Weibo He.
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

#include <Physica/Core/Exception/CUDA/CUDA.cuh>

namespace Physica::Core {
    class PHYSICA_API CudaEvent {
        cudaEvent_t event;
    public:
        CudaEvent();
        CudaEvent(const CudaEvent&) = delete;
        CudaEvent(CudaEvent&& obj) noexcept;
        ~CudaEvent();
        /* Operators */
        CudaEvent& operator=(CudaEvent obj) noexcept;
        /* Operations */
        inline void wait();
        void swap(CudaEvent& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] cudaEvent_t getEvent() const noexcept { return event; }
    };

    inline void CudaEvent::wait() {
        check(cudaEventSynchronize(event));
    }
}
