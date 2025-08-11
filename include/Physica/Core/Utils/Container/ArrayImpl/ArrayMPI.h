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
#pragma once

#include "../Array.h"
#include "Physica/Core/Parallel/MPIContext.h"
#include "Physica/Core/Exception/MPIException.h"

namespace Physica {
    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::send(int from, int to) {
        assert(from < MPIContext::getNumProcess() && to < MPIContext::getNumProcess() && from != to);
        const int id = MPIContext::getProcessID();
        if (id == from)
            check_mpi(MPI_Send(data(), getLength(), T::dtype_mpi(), to, 0, MPIContext::getWorld()));
        else if (id == to)
            check_mpi(MPI_Recv(data(), getLength(), T::dtype_mpi(), from, MPI_ANY_TAG, MPIContext::getWorld(), MPI_STATUS_IGNORE));
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::sendrecv(int send_to, int recv_from) {
        check_mpi(MPI_Sendrecv_replace(
                data(), getLength(), T::dtype_mpi(), send_to, 0, recv_from, MPI_ANY_TAG, MPIContext::getWorld(), MPI_STATUS_IGNORE));
    }
}
