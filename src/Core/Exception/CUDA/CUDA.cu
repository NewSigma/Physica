/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <system_error>
#include "../SystemError.h"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"

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
        [[nodiscard]] const char* name() const noexcept final { return "CUDA Runtime"; }
        [[nodiscard]] std::string message(int code) const final { return cudaGetErrorString(cudaError_t(code)); }
    };
}

void Physica::throw_cuda_error(int code) {
    constinit static const Impl instance{};
    throw SystemError(code, instance);
}
