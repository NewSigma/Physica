/*
 * Copyright 2022 WeiBo He.
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

#include <exception>
#include <cuda_runtime.h>
#include "Physica/Macro.h"

namespace Physica::Core {
    class PHYSICA_API CudaException : public std::exception {
        cudaError_t code;
    public:
        CudaException(cudaError_t code_) : code(code_) {}
        ~CudaException() noexcept override = default;
        /* Getters */
        [[nodiscard]] const char* what() const noexcept override { return cudaGetErrorString(code); }
        [[nodiscard]] cudaError_t getCode() const noexcept { return code; }
    };
}
