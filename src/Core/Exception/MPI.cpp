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
#include <array>
#include <format>
#include <system_error>
#include "Physica/Core/Exception/MPI.h"
#include "Physica/Core/Parallel/MPI.h"
#include <mpi.h>

using namespace Physica;

namespace {
    class Impl final : public std::error_category {
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
}

std::string Impl::message(int err) const {
    std::array<char, MPI_MAX_ERROR_STRING> buffer{};
    int resultlen{};
    MPI_Error_string(err, buffer.data(), &resultlen);
    return std::format("Rank {}: {}", MPI::getRank(), buffer.data());
}

void Physica::check_mpi(int err) {
    constinit static const Impl instance{};
    if (err != 0) [[unlikely]]
        throw std::system_error(err, instance);
}
