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
#include <iostream>
#include <format>
#include "Physica/Core/Parallel/MPI.h"
#include "Physica/Core/Exception/MPIException.h"
#include "Physica/Core/Utils/Builtin.h"
#include "Physica/Core/Utils/NoImpl.h"
#ifdef PHYSICA_MPI
    #include <mpi/mpi.h>
#else
constexpr static void* MPI_COMM_WORLD = nullptr;
#endif

using namespace Physica;

namespace {
    constexpr void checkPID([[maybe_unused]] int pid) noexcept {
        assert(0 <= pid && pid < MPI::getNumProcess());
    }

    constexpr void checkPID(int pid, auto... pids) {
        checkPID(pid);
        checkPID(pids...);
    }
}

const auto MPI::World = Handle<HandleType::MPI_Comm>(MPI_COMM_WORLD);

MPI::MPI() noexcept {
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

MPI::~MPI() {
#ifdef PHYSICA_MPI
    MPI_Finalize();
#endif
}

MPI& MPI::getInstance() noexcept {
    static MPI mpi{};
    return mpi;
}

bool MPI::initialized() noexcept {
#ifdef PHYSICA_MPI
    int result = 0;
    MPI_Initialized(&result);
    return result != 0;
#else
    return false;
#endif
}

int MPI::getNumProcess() noexcept {
#ifdef PHYSICA_MPI
    int result{};
    MPI_Comm_size(MPI_Comm(World), &result);
    return result;
#else
    return 1;
#endif
}

int MPI::getProcessID() noexcept {
#ifdef PHYSICA_MPI
    int result{};
    MPI_Comm_rank(MPI_Comm(World), &result);
    return result;
#else
    return 0;
#endif
}

void MPI::send(int to, void* data, int count, dtype_handle dtype, comm_handle comm) {
    checkPID(to);
#ifdef PHYSICA_MPI

    check_mpi(MPI_Send(data, count, MPI_Datatype(dtype), to, 0, MPI_Comm(comm)));
#endif
}

void MPI::recv(int from, void* data, int count, dtype_handle dtype, comm_handle comm) {
    checkPID(from);
#ifdef PHYSICA_MPI
    MPI_Status status;
    check_mpi(MPI_Recv(data, count, MPI_Datatype(dtype), from, MPI_ANY_TAG, MPI_Comm(comm), IsDebug() ? &status : MPI_STATUS_IGNORE));
    if constexpr (IsDebug()) {
        int recv_count{};
        MPI_Get_count(&status, MPI_Datatype(dtype), &recv_count);
        assert(recv_count == count && "[Error]: Recv length do not match");
    }
#endif
}
/**
 * A paired send-recv, useful for avoiding dead lock
 */
void MPI::pass(int from, int to, void* data, int count, dtype_handle dtype, comm_handle comm) {
    assert(from != to);
    const int id = MPI::getProcessID();
    if (id == from)
        send(to, data, count, dtype, comm);
    else if (id == to)
        recv(from, data, count, dtype, comm);
}

void MPI::sendrecv(int send_to, int recv_from, void* data, int count, dtype_handle dtype, comm_handle comm) {
    checkPID(send_to, recv_from);
#ifdef PHYSICA_MPI
    MPI_Status status;
    check_mpi(MPI_Sendrecv_replace(data, count, MPI_Datatype(dtype), send_to, 0, recv_from, MPI_ANY_TAG, MPI_Comm(comm), IsDebug() ? &status : MPI_STATUS_IGNORE));
    if constexpr (IsDebug()) {
        int recv_count{};
        MPI_Get_count(&status, MPI_Datatype(dtype), &recv_count);
        assert(recv_count == count && "[Error]: Recv length do not match");
    }
#endif
}

void MPI::bcast(int root, void* data, int count, dtype_handle dtype, comm_handle comm) {
    checkPID(root);
#ifdef PHYSICA_MPI
    check_mpi(MPI_Bcast(data, count, MPI_Datatype(dtype), root, MPI_Comm(comm)));
#endif
}

void MPI::wait(comm_handle comm) {
#ifdef PHYSICA_MPI
    check_mpi(MPI_Barrier(MPI_Comm(comm)));
#endif
}

auto MPI::dtype_primitive(PrimitiveType type) noexcept -> dtype_handle {
    switch (type) {
#ifdef PHYSICA_MPI
    case PrimitiveType::Int8:
        return dtype_handle(MPI_INT8_T);
    case PrimitiveType::Int16:
        return dtype_handle(MPI_INT16_T);
    case PrimitiveType::Int32:
        return dtype_handle(MPI_INT32_T);
    case PrimitiveType::Int64:
        return dtype_handle(MPI_INT64_T);
    case PrimitiveType::UInt8:
        return dtype_handle(MPI_UINT8_T);
    case PrimitiveType::UInt16:
        return dtype_handle(MPI_UINT16_T);
    case PrimitiveType::UInt32:
        return dtype_handle(MPI_UINT32_T);
    case PrimitiveType::UInt64:
        return dtype_handle(MPI_UINT64_T);
    case PrimitiveType::Bool:
        return dtype_handle(MPI_C_BOOL);
    case PrimitiveType::Char:
        return dtype_handle(MPI_CHAR);
    case PrimitiveType::SignedChar:
        return dtype_handle(MPI_SIGNED_CHAR);
    case PrimitiveType::UnsignedChar:
        return dtype_handle(MPI_UNSIGNED_CHAR);
    case PrimitiveType::Short:
        return dtype_handle(MPI_SHORT);
    case PrimitiveType::UnsignedShort:
        return dtype_handle(MPI_UNSIGNED_SHORT);
    case PrimitiveType::Int:
        return dtype_handle(MPI_INT);
    case PrimitiveType::UnsignedInt:
        return dtype_handle(MPI_UNSIGNED);
    case PrimitiveType::Long:
        return dtype_handle(MPI_LONG);
    case PrimitiveType::UnsignedLong:
        return dtype_handle(MPI_UNSIGNED_LONG);
    case PrimitiveType::LongLong:
        return dtype_handle(MPI_LONG_LONG);
    case PrimitiveType::UnsignedLongLong:
        return dtype_handle(MPI_UNSIGNED_LONG_LONG);
    case PrimitiveType::Float:
        return dtype_handle(MPI_FLOAT);
    case PrimitiveType::Double:
        return dtype_handle(MPI_DOUBLE);
    case PrimitiveType::LongDouble:
        return dtype_handle(MPI_LONG_DOUBLE);
#endif
    default:
        unreachable();
    }
}
