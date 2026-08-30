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
#include <format>
#include "Physica/Core/Parallel/MPIContext.h"
#include "Physica/Core/Exception/MPIException.h"
#include "Physica/Core/Utils/NoImpl.h"
#ifdef PHYSICA_MPI
    #include <mpi/mpi.h>
#else
    constexpr void* MPI_COMM_WORLD = nullptr;
#endif

using namespace Physica;

const auto MPIContext::World = Handle<HandleType::MPI_Comm>(MPI_COMM_WORLD);

MPIContext::MPIContext() noexcept {
#ifdef PHYSICA_MPI
    try {
        int mode = -1;
        check_mpi(MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED, &mode));
        if (mode != MPI_THREAD_SERIALIZED) {
            std::cerr << "[Error]: Physica do not support the MPI\n";
            std::abort();
        }
        check_mpi(MPI_Comm_set_errhandler(MPI_Comm(World), MPI_ERRORS_RETURN));
    }
    catch (std::exception& e) {
        std::cerr << std::format("MPI init failed: {}\n", e.what());
        MPI_Abort(MPI_Comm(World), -1);
    }
#else
    noImpl("[Error]: Cannot initialize MPI if not compiled with MPI support");
#endif
}

MPIContext::~MPIContext() {
#ifdef PHYSICA_MPI
    MPI_Finalize();
#endif
}

MPIContext& MPIContext::getInstance() noexcept {
    static MPIContext mpi{};
    return mpi;
}

int MPIContext::getNumProcess() noexcept {
#ifdef PHYSICA_MPI
    int result{};
    MPI_Comm_size(MPI_Comm(World), &result);
    return result;
#else
    return 1;
#endif
}

int MPIContext::getProcessID() noexcept {
#ifdef PHYSICA_MPI
    int result{};
    MPI_Comm_rank(MPI_Comm(World), &result);
    return result;
#else
    return 0;
#endif
}

void MPIContext::send(int from, int to, void* data, int count, dtype_handle dtype, comm_handle comm) {
#ifdef PHYSICA_MPI
    assert(0 < from && from < MPIContext::getNumProcess());
    assert(0 < to && to < MPIContext::getNumProcess());
    assert(from != to);
    const int id = MPIContext::getProcessID();
    if (id == from)
        check_mpi(MPI_Send(data, count, MPI_Datatype(dtype), to, 0, MPI_Comm(comm)));
    else if (id == to) {
        MPI_Status status;
        check_mpi(MPI_Recv(data, count, MPI_Datatype(dtype), from, MPI_ANY_TAG, MPI_Comm(comm), IsDebug() ? &status : MPI_STATUS_IGNORE));
        if constexpr (IsDebug()) {
            int recv_count{};
            MPI_Get_count(&status, MPI_Datatype(dtype), &recv_count);
            assert(recv_count == count && "[Error]: Recv length do not match");
        }
    }
#endif
}

void MPIContext::sendrecv(int send_to, int recv_from, void* data, int count, dtype_handle dtype, comm_handle comm) {
#ifdef PHYSICA_MPI
    assert(0 < send_to && send_to < MPIContext::getNumProcess());
    assert(0 < recv_from && recv_from < MPIContext::getNumProcess());
    MPI_Status status;
    check_mpi(MPI_Sendrecv_replace(data, count, MPI_Datatype(dtype), send_to, 0, recv_from, MPI_ANY_TAG, MPI_Comm(comm), IsDebug() ? &status : MPI_STATUS_IGNORE));
    if constexpr (IsDebug()) {
        int recv_count{};
        MPI_Get_count(&status, MPI_Datatype(dtype), &recv_count);
        assert(recv_count == count && "[Error]: Recv length do not match");
    }
#endif
}

void MPIContext::bcast(int root, void* data, int count, dtype_handle dtype, comm_handle comm) {
#ifdef PHYSICA_MPI
    assert(0 < root && root < MPIContext::getNumProcess());
    check_mpi(MPI_Bcast(data, count, MPI_Datatype(dtype), root, MPI_Comm(comm)));
#endif
}

void MPIContext::wait(comm_handle comm) {
#ifdef PHYSICA_MPI
    check_mpi(MPI_Barrier(MPI_Comm(comm)));
#endif
}
