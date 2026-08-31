/*
 * Copyright 2025-2026 Weibo He.
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
#include "Test.h"
#include <cstdlib>
#include <print>
#include "Physica/Core/Parallel/MPIContext.h"
#ifdef PHYSICA_MPI
    #include <mpi/mpi.h>
#endif

namespace {
    [[noreturn, gnu::cold]] void expect_fail(std::source_location loc, uint64_t seed) noexcept {
        std::println("Failed at file: {}:{}:{}", loc.file_name(), loc.line(), loc.column());
        std::println("          func: {}", loc.function_name());
        if (seed != 0)
            std::println("          seed: {}", seed);
    #ifdef PHYSICA_MPI
        using namespace Physica;
        if (MPIContext::initialized())
            MPI_Abort(MPI_Comm(MPIContext::World), EXIT_FAILURE);
    #endif
        exit(EXIT_FAILURE);
    }
}

namespace Physica {
    void expect(bool pass, std::source_location loc, uint64_t seed) noexcept {
        if (!pass) [[unlikely]]
            expect_fail(loc, seed);
    }
}
