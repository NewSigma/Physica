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
#include "Physica/Core/Exception/MPIException.h"
#include "Physica/Core/Parallel/MPI.h"
#include <array>
#include <format>
#include <mpi/mpi.h>

using namespace Physica;

namespace {
    class Impl final : public std::error_category {
        int rank = MPI::getRank();
    public:
        Impl() = default;
        Impl(const Impl&) = delete;
        Impl(Impl&&) noexcept = delete;
        ~Impl() = default;
        /* Operators */
        Impl& operator=(const Impl&) = delete;
        Impl& operator=(Impl&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] const char* name() const noexcept final { return "MPI"; }
        [[nodiscard]] std::string message(int) const final;
    };

    std::string Impl::message(int err) const {
        std::array<char, MPI_MAX_ERROR_STRING> buffer{};
        int resultlen{};
        MPI_Error_string(err, buffer.data(), &resultlen);
        return std::format("Rank {}: {}", rank, buffer.data());
    }
}

MPIException::MPIException(int code) noexcept : std::system_error(code, Impl()) {}
