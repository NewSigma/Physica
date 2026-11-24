/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Test.h"
#include <cstdlib>
#include <print>
#include "Physica/Core/Parallel/MPI.h"
#ifdef PHYSICA_MPI
    #include <mpi.h>
#endif

namespace {
    [[noreturn, gnu::cold]] void expect_fail(std::source_location loc, uint64_t seed) noexcept {
        std::println("Failed at file: {}:{}:{}", loc.file_name(), loc.line(), loc.column());
        std::println("          func: {}", loc.function_name());
        if (seed != 0)
            std::println("          seed: {}", seed);
        std::fflush(nullptr);

    #ifdef PHYSICA_MPI
        using namespace Physica;
        if (MPI::initialized() && MPI::getRank() > 1)
            MPI_Abort(MPI_Comm(MPI::World), EXIT_FAILURE);
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
