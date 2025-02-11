/*
 * Copyright 2024 Weibo He.
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
#include <format>
#include <sstream>
#ifdef PHYSICA_MKL
    #include <mkl_service.h>
#endif
#ifdef PHYSICA_CUDA
    #include "Physica/Core/Exception/CUDA/CUDA.cuh"
#endif
#include "Physica/Macro.h"

namespace Physica {
    PHYSICA_API std::string version() {
        std::ostringstream os{};
        std::format_to(std::ostreambuf_iterator<char>(os), "Physica {}    Hash: {}\n", Version, GitHash);
    #ifdef PHYSICA_MKL
        constexpr int len = 198;
        char buf[198];
        mkl_get_version_string(buf, len);
        std::format_to(std::ostreambuf_iterator<char>(os), "MKL: {}\n", buf);
    #endif
    #ifdef PHYSICA_CUDA
        int driverVer, runtimeVer;
        check(cudaDriverGetVersion(&driverVer));
        check(cudaRuntimeGetVersion(&runtimeVer));
        std::format_to(std::ostreambuf_iterator<char>(os), "CUDA driver : {}\n", driverVer);
        std::format_to(std::ostreambuf_iterator<char>(os), "CUDA runtime: {}\n", runtimeVer);
    #endif
        return os.str();
    }
}
