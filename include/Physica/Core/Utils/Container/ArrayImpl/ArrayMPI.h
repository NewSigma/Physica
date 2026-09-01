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
#include "Physica/Core/Parallel/MPI.h"
#include "Physica/Core/Exception/MPIException.h"

namespace Physica {
    template<class Derived, class Allocator>
    void ArrayBase<Derived, Allocator>::pass(int from, int to) {
        MPI::pass(from, to, data(), getLength(), MPI::dtype<ElemType>(), MPI::World);
    }

    template<class Derived, class Allocator>
    void ArrayBase<Derived, Allocator>::sendrecv(int send_to, int recv_from) {
        MPI::sendrecv(send_to, recv_from, data(), getLength(), MPI::dtype<ElemType>(), MPI::World);
    }

    template<class Derived, class Allocator>
    void ArrayBase<Derived, Allocator>::bcast(int root) {
        MPI::bcast(root, data(), getLength(), MPI::dtype<ElemType>(), MPI::World);
    }
}
