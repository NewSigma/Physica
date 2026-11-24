/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <format>
#include "Physica/Macro.h"

namespace Physica {
    class KernelConfig {
        using This = KernelConfig;
    public:
        dim3 blocks;
        dim3 threads;
    public:
        __host__ __device__ constexpr KernelConfig(dim3 blocks_, dim3 threads_) noexcept : blocks(blocks_), threads(threads_) {}
        constexpr KernelConfig(const This&) = default;
        constexpr KernelConfig(This&&) noexcept = default;
        constexpr ~KernelConfig() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        PHYSICA_API void dump();
    };
}

namespace std {
    template<>
    struct PHYSICA_API formatter<Physica::KernelConfig, char> {
        constexpr static auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
        static auto format(const Physica::KernelConfig& obj, std::format_context& ctx) -> std::format_context::iterator;
    };
}
