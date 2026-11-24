/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <concepts>

namespace Physica {
    enum ExecutePolicy : int8_t {
        Sequential,
        Thread,
        GPU,
    };

    class EmptyTask;
    class Task;

    template<ExecutePolicy P>
    [[nodiscard]] auto schedule(std::invocable<> auto fn);

    template<ExecutePolicy P>
    [[nodiscard]] auto parallel_for(auto fn, size_t num_loop);

    template<ExecutePolicy P>
    [[nodiscard]] auto parallel_for(auto fn, size_t num_loop, size_t part);
}
