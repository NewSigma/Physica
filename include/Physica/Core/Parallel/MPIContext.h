/*
 * Copyright 2024-2026 Weibo He.
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

#ifdef PHYSICA_MPI
    #include <mpi/mpi.h>
#endif
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API MPIContext final {
        using This = MPIContext;
    public:
        MPIContext(const This&) = delete;
        MPIContext(This&&) noexcept = delete;
        ~MPIContext();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Static memebers */
        [[nodiscard]] static MPIContext& getInstance() noexcept;
        [[nodiscard]] constexpr static auto getWorld() noexcept;
        [[nodiscard]] static inline int getNumProcess() noexcept;
        [[nodiscard]] static inline int getProcessID() noexcept;
        static void wait();
    private:
        MPIContext() noexcept;
    };

    constexpr auto MPIContext::getWorld() noexcept {
    #ifdef PHYSICA_MPI
        return MPI_COMM_WORLD;
    #else
        return 0;
    #endif
    }

    int MPIContext::getNumProcess() noexcept {
        int result = 1;
    #ifdef PHYSICA_MPI
        MPI_Comm_size(getWorld(), &result);
    #endif
        return result;
    }

    int MPIContext::getProcessID() noexcept {
        int result = 0;
    #ifdef PHYSICA_MPI
        MPI_Comm_rank(getWorld(), &result);
    #endif
        return result;
    }
}
