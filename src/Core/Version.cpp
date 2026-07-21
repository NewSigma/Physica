/*
 * Copyright 2024-2025 Weibo He.
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
module;

#include <format>
#include <sstream>
#ifdef PHYSICA_HDF5
    #include <hdf5.h>
#endif
#ifdef PHYSICA_MKL
    #include <mkl_service.h>
#endif
#ifdef PHYSICA_MPI
    #include <mpi.h>
#endif
#ifdef PHYSICA_CUDA
    #include "Physica/Core/Exception/CUDA/CUDA.cuh"
#endif
#ifdef PHYSICA_MIMALLOC
    #include <mimalloc.h>
#endif
#include "Physica/Macro.h"

module Physica.Core;

PHYSICA_API std::string Physica::version() {
    std::ostringstream os{};
    std::format_to(std::ostreambuf_iterator<char>(os), "Physica {}    Hash: {}\n", Version, GitHash);
    {
    #ifdef PHYSICA_HDF5
        unsigned int major{};
        unsigned int minor{};
        unsigned int release{};
        H5get_libversion(&major, &minor, &release);
        std::format_to(std::ostreambuf_iterator<char>(os), "HDF5: {}.{}-{}\n", major, minor, release);
    #endif
    }
    {
    #ifdef PHYSICA_MKL
        std::array<char, 198> buf{}; // 198 referenced from MKL Doc
        mkl_get_version_string(buf.data(), buf.size());
        std::format_to(std::ostreambuf_iterator<char>(os), "MKL: {}\n", buf.data());
    #endif
    }
    {
    #ifdef PHYSICA_MPI
        std::array<char, MPI_MAX_LIBRARY_VERSION_STRING> buf{};
        int len{};
        MPI_Get_library_version(buf.data(), &len);
        std::format_to(std::ostreambuf_iterator<char>(os), "MPI: {}\n", buf.data());
    #endif
    }
    {
    #ifdef PHYSICA_CUDA
        int driverVer{};
        int runtimeVer{};
        check(cudaDriverGetVersion(&driverVer));
        check(cudaRuntimeGetVersion(&runtimeVer));
        std::format_to(std::ostreambuf_iterator<char>(os), "CUDA driver : {}\n", driverVer);
        std::format_to(std::ostreambuf_iterator<char>(os), "CUDA runtime: {}\n", runtimeVer);
    #endif
    }
#ifdef PHYSICA_MIMALLOC
    std::format_to(std::ostreambuf_iterator<char>(os), "mimalloc: {}\n", mi_version());
#endif
    return os.str();
}
