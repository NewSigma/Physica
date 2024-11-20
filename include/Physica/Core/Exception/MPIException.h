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

#include <exception>
#include <mpi/mpi.h>
#include "Physica/Macro.h"

namespace Physica::Core {
    class PHYSICA_API MPIException : public std::exception {
        char* msg;
    public:
        MPIException(const char* msg_) noexcept;
        MPIException(int err) noexcept;
        ~MPIException() override;
        const char* what() const noexcept override { return msg; }
    };
}

namespace Physica {
    inline void mpiCheck(int err) {
        if (err != MPI_SUCCESS) [[unlikely]]
            throw Physica::Core::MPIException(err);
    }
}

#endif
