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
#include <cassert>
#include "Physica/Core/Parallel/CUDAEvent.cuh"

namespace Physica::Core {
    CUDAEvent::CUDAEvent() {
        check(cudaEventCreate(&event));
    }

    CUDAEvent::CUDAEvent(CUDAEvent&& obj) noexcept : event(obj.event) {
        obj.event = nullptr;
    }

    CUDAEvent::~CUDAEvent() {
        cudaEventDestroy(event);
    }

    void CUDAEvent::swap(CUDAEvent& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(event, obj.event);
    }
}
