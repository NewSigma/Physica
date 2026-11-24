/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <print>
#include "Physica/Core/Parallel/KernelConfig.cuh"

using namespace Physica;

void KernelConfig::dump() {
    std::println("{}", *this);
}

namespace std {
    auto formatter<KernelConfig, char>::format(
            const KernelConfig& obj, std::format_context& ctx) -> std::format_context::iterator{
        return std::format_to(ctx.out(), "Blocks: ({}, {}, {})\nThreads: ({}, {}, {})",
                obj.blocks.x, obj.blocks.y, obj.blocks.z, obj.threads.x, obj.threads.y, obj.threads.z);
    }
}
