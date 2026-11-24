/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <cassert>
#include <utility>
#include "Physica/Core/Parallel/MPI/Request.h"
#include "Physica/Core/Exception/MPI.h"
#ifdef PHYSICA_MPI
    #include <mpi.h>
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
