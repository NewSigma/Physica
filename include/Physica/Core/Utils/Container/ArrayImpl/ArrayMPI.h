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

#include "ArrayBase.h"
#include "Physica/Core/Parallel/MPIContext.h"
#include "Physica/Core/Exception/MPIException.h"

namespace Physica {
    template<class Derived, class Allocator>
    void ArrayBase<Derived, Allocator>::send(int from, int to) {
        assert(0 < from && from < MPIContext::getNumProcess());
        assert(0 < to && to < MPIContext::getNumProcess());
        assert(from != to);
        const int id = MPIContext::getProcessID();
        if (id == from)
            check_mpi(MPI_Send(data(), getLength(), ElemType::dtype_mpi(), to, 0, MPIContext::getWorld()));
        else if (id == to) {
            MPI_Status status;
            check_mpi(MPI_Recv(data(), getLength(), ElemType::dtype_mpi(), from, MPI_ANY_TAG, MPIContext::getWorld(), IsDebug() ? &status : MPI_STATUS_IGNORE));
            if constexpr (IsDebug()) {
                int count{};
                MPI_Get_count(&status, ElemType::dtype_mpi(), &count);
                assert(count == getLength() && "[Error]: Recv length do not match");
            }
        }
    }

    template<class Derived, class Allocator>
    void ArrayBase<Derived, Allocator>::sendrecv(int send_to, int recv_from) {
        assert(0 < send_to && send_to < MPIContext::getNumProcess());
        assert(0 < recv_from && recv_from < MPIContext::getNumProcess());
        MPI_Status status;
        check_mpi(MPI_Sendrecv_replace(data(), getLength(), ElemType::dtype_mpi(), send_to, 0, recv_from, MPI_ANY_TAG, MPIContext::getWorld(), IsDebug() ? &status : MPI_STATUS_IGNORE));
        if constexpr (IsDebug()) {
            int count{};
            MPI_Get_count(&status, ElemType::dtype_mpi(), &count);
            assert(count == getLength() && "[Error]: Recv length do not match");
        }
    }

    template<class Derived, class Allocator>
    void ArrayBase<Derived, Allocator>::bcast(int root) {
        assert(0 < root && root < MPIContext::getNumProcess());
        check_mpi(MPI_Bcast(data(), getLength(), ElemType::dtype_mpi(), root, MPIContext::getWorld()));
    }
}
