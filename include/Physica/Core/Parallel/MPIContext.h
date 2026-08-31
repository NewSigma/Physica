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
#pragma once

#include "Physica/Macro.h"
#include "Physica/Core/Utils/Handle.h"

namespace Physica {
    class PHYSICA_API MPIContext final {
        using This = MPIContext;
    public:
        using comm_handle = Handle<HandleType::MPI_Comm>;
        using dtype_handle = Handle<HandleType::MPI_Dtype>;

        static const comm_handle World;
    public:
        MPIContext(const This&) = delete;
        MPIContext(This&&) noexcept = delete;
        ~MPIContext();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Static memebers */
        [[nodiscard]] static This& getInstance() noexcept;
        [[nodiscard]] static bool initialized() noexcept;
        [[nodiscard]] static int getNumProcess() noexcept;
        [[nodiscard]] static int getProcessID() noexcept;

        static void send(int from, int to, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void sendrecv(int send_to, int recv_from, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void bcast(int root, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void wait(comm_handle comm = World);
    private:
        MPIContext() noexcept;
    };
}
