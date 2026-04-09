/*
 * Copyright 2025-2026 Weibo He.
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
#include <print>
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"

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
