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

#include <format>
#include "Physica/Core/Parallel/MPIContext.h"
#include "Physica/Core/Exception/MPIException.h"

using namespace Physica;

MPIContext::MPIContext() noexcept {
    try {
        int mode = -1;
        check_mpi(MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED, &mode));
        if (mode != MPI_THREAD_SERIALIZED) {
            std::cerr << "[Error]: Physica do not support the MPI\n";
            std::abort();
        }
        check_mpi(MPI_Comm_set_errhandler(getWorld(), MPI_ERRORS_RETURN));
    }
    catch (std::exception& e) {
        std::cerr << std::format("MPI init failed: {}\n", e.what());
        MPI_Abort(getWorld(), -1);
    }
}

MPIContext::~MPIContext() {
    MPI_Finalize();
}

MPIContext& MPIContext::getInstance() noexcept {
    static MPIContext mpi{};
    return mpi;
}

void MPIContext::wait() {
    check_mpi(MPI_Barrier(getWorld()));
}

#endif
