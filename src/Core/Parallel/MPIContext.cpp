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
#include "Physica/Core/Utils/Builtin.h"
#include "Physica/Core/Utils/NoImpl.h"
#ifdef PHYSICA_MPI
    #include <mpi/mpi.h>
#else
    constexpr void* MPI_COMM_WORLD = nullptr;
#endif

using namespace Physica;

namespace {
    [[maybe_unused]] constexpr void checkPID([[maybe_unused]] int pid) noexcept {
        assert(0 <= pid && pid < MPIContext::getNumProcess());
    }

    [[maybe_unused]] constexpr void checkPID(int pid, auto... pids) {
        checkPID(pid);
        checkPID(pids...);
    }
}

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
        MPI_Abort(MPI_Comm(World), EXIT_FAILURE);
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

bool MPIContext::initialized() noexcept {
#ifdef PHYSICA_MPI
    int result = 0;
    MPI_Initialized(&result);
    return result != 0;
#else
    return false;
#endif
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
    checkPID(from, to);
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
    checkPID(send_to, recv_from);
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
    checkPID(root);
    check_mpi(MPI_Bcast(data, count, MPI_Datatype(dtype), root, MPI_Comm(comm)));
#endif
}

void MPIContext::wait(comm_handle comm) {
#ifdef PHYSICA_MPI
    check_mpi(MPI_Barrier(MPI_Comm(comm)));
#endif
}

auto MPIContext::dtype_scalar(Dtype type) noexcept -> dtype_handle {
    switch (type) {
#ifdef PHYSICA_MPI
    case Dtype::Int8:
        return dtype_handle(MPI_INT8_T);
    case Dtype::Int16:
        return dtype_handle(MPI_INT16_T);
    case Dtype::Int32:
        return dtype_handle(MPI_INT32_T);
    case Dtype::Int64:
        return dtype_handle(MPI_INT64_T);
    case Dtype::UInt8:
        return dtype_handle(MPI_UINT8_T);
    case Dtype::UInt16:
        return dtype_handle(MPI_UINT16_T);
    case Dtype::UInt32:
        return dtype_handle(MPI_UINT32_T);
    case Dtype::UInt64:
        return dtype_handle(MPI_UINT64_T);
    case Dtype::Bool:
        return dtype_handle(MPI_C_BOOL);
    case Dtype::Char:
        return dtype_handle(MPI_CHAR);
    case Dtype::SignedChar:
        return dtype_handle(MPI_SIGNED_CHAR);
    case Dtype::UnsignedChar:
        return dtype_handle(MPI_UNSIGNED_CHAR);
    case Dtype::Short:
        return dtype_handle(MPI_SHORT);
    case Dtype::UnsignedShort:
        return dtype_handle(MPI_UNSIGNED_SHORT);
    case Dtype::Int:
        return dtype_handle(MPI_INT);
    case Dtype::UnsignedInt:
        return dtype_handle(MPI_UNSIGNED);
    case Dtype::Long:
        return dtype_handle(MPI_LONG);
    case Dtype::UnsignedLong:
        return dtype_handle(MPI_UNSIGNED_LONG);
    case Dtype::LongLong:
        return dtype_handle(MPI_LONG_LONG);
    case Dtype::UnsignedLongLong:
        return dtype_handle(MPI_UNSIGNED_LONG_LONG);
    case Dtype::Float:
        return dtype_handle(MPI_FLOAT);
    case Dtype::Double:
        return dtype_handle(MPI_DOUBLE);
    case Dtype::LongDouble:
        return dtype_handle(MPI_LONG_DOUBLE);
#endif
    default:
        unreachable();
    }
}
