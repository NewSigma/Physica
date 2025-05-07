/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Exception/CUDA/cuRAND.cuh"

using namespace Physica;

std::string cuRANDException::Impl::message(int code) const {
    switch (code) {
    case CURAND_STATUS_SUCCESS:;
        return "No error";
    case CURAND_STATUS_VERSION_MISMATCH:
        return "Header file and linked library version do not match";
    case CURAND_STATUS_NOT_INITIALIZED:
        return "Generator not initialized";
    case CURAND_STATUS_ALLOCATION_FAILED:
        return "Memory allocation failed";
    case CURAND_STATUS_TYPE_ERROR:
        return "Generator is wrong type";
    case CURAND_STATUS_OUT_OF_RANGE:
        return "Argument out of range";
    case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
        return "Length requested is not a multple of dimension";
    case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
        return "GPU does not have double precision required by MRG32k3a";
    case CURAND_STATUS_LAUNCH_FAILURE:
        return "Kernel launch failure";
    case CURAND_STATUS_PREEXISTING_FAILURE:
        return "Preexisting failure on library entry";
    case CURAND_STATUS_INITIALIZATION_FAILED:
        return "Initialization of CUDA failed";
    case CURAND_STATUS_ARCH_MISMATCH:
        return "Architecture mismatch, GPU does not support requested feature";
    case CURAND_STATUS_INTERNAL_ERROR:
        return "Internal library error";
    default:
        return "Bad error code";
    }
}
