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
