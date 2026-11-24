/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <array>
#include <format>
#include <system_error>
#include "SystemError.h"
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
        throw SystemError(err, instance);
}
