/*
 * Copyright 2026 Weibo He.
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
#include <cassert>
#include <utility>
#include "Physica/Core/Parallel/MPI/Request.h"
#include "Physica/Core/Exception/MPIException.h"
#ifdef PHYSICA_MPI
    #include <mpi/mpi.h>
#endif

using namespace Physica;

namespace {
    auto empty_request() {
    #ifdef PHYSICA_MPI
        return MPI_REQUEST_NULL;
    #else
        return nullptr;
    #endif
    }
}

Request::Request() : h(empty_request()) {}

Request::Request(Handle h) : h(h) {}

Request::Request(Request&& obj) noexcept : h(std::exchange(obj.h, Handle(empty_request()))) {}

Request::~Request() {
#ifdef PHYSICA_MPI
    MPI_Request req{};
    req = MPI_Request(h);
    if (req != empty_request())
        MPI_Request_free(&req);
#endif
}

void Request::wait() {
#ifdef PHYSICA_MPI
    MPI_Request req{};
    req = MPI_Request(h);
    check_mpi(MPI_Wait(&req, MPI_STATUS_IGNORE));
    h = Handle(req);
#endif
}

void Request::swap(Request& obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    h.swap(obj.h);
}

bool Request::query() {
#ifdef PHYSICA_MPI
    MPI_Request req{};
    req = MPI_Request(h);
    int flag = 0;
    check_mpi(MPI_Test(&req, &flag, MPI_STATUS_IGNORE));
    h = Handle(req);
    return flag != 0;
#else
    return false;
#endif
}
