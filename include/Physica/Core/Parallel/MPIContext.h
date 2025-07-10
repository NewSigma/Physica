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
#pragma once

#ifdef PHYSICA_MPI

#include <mpi/mpi.h>
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API MPIContext {
        using This = MPIContext;
    public:
        MPIContext(const This&) = delete;
        MPIContext(This&&) noexcept = delete;
        ~MPIContext();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] inline int getNumProcess() const noexcept;
        [[nodiscard]] inline int getProcessID() const noexcept; 
        /* Static memebers */
        [[nodiscard]] static MPIContext& getInstance() noexcept;
    private:
        MPIContext();
        /* Getters */
        static MPI_Comm getWorld() noexcept { return MPI_COMM_WORLD; }
        /* Static memebers */
        static void world_handler(MPI_Comm* pComm, int* pErr, ...);
    };

    inline int MPIContext::getNumProcess() const noexcept {
        int result;
        MPI_Comm_size(getWorld(), &result);
        return result;
    }

    inline int MPIContext::getProcessID() const noexcept {
        int result;
        MPI_Comm_rank(getWorld(), &result);
        return result;
    }
}

#endif
