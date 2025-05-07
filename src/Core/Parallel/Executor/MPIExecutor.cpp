/*
 * Copyright 2024-2025 Weibo He.
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
#ifdef PHYSICA_MPI

#include "Physica/Core/Parallel/Executor/MPIExecutor.h"
#include "Physica/Core/Exception/MPIException.h"

using namespace Physica;

MPIExecutor::MPIExecutor() {
    int mode;
    check_mpi(MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED, &mode));
    if (mode != MPI_THREAD_SERIALIZED)
        throw std::runtime_error("[Error]: Physica do not support the MPI");

    MPI_Errhandler handler;
    check_mpi(MPI_Comm_create_errhandler(world_handler, &handler));
    check_mpi(MPI_Comm_set_errhandler(getWorld(), handler));
}

MPIExecutor::~MPIExecutor() {
    MPI_Finalize();
}

MPIExecutor& MPIExecutor::getInstance() noexcept {
    static MPIExecutor mpi{};
    return mpi;
}

void MPIExecutor::world_handler(MPI_Comm*, int* pErr, ...) {
    throw MPIException(*pErr);
}

#endif
