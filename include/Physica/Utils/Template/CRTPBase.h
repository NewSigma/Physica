/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Config.h"
#ifdef PHYSICA_CUDA
    #include <cuda_runtime_api.h>
#endif

namespace Physica::Utils {
    /**
     * This class helps implemant CRTP.
     */
    template<class Derived>
    class CRTPBase {
    public:
        [[nodiscard]] __host__ __device__ Derived& getDerived() noexcept { return *static_cast<Derived*>(this); }
        [[nodiscard]] __host__ __device__ const Derived& getDerived() const noexcept { return *static_cast<const Derived*>(this); }
        [[nodiscard]] __host__ __device__ const Derived& getConstDerived() const noexcept { return *static_cast<Derived*>(this); }
        [[nodiscard]] __host__ __device__ Derived& getConstCastDerived() const noexcept { return *static_cast<Derived*>(const_cast<CRTPBase*>(this)); }
    };
}
