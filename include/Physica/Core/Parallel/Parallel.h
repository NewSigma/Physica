/*
 * Copyright 2022-2025 Weibo He.
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

#include "Task.h"

namespace Physica {
    template<ParallelPolicy P>
    [[nodiscard]] Task<P> schedule(auto func) noexcept {
        func();
        co_return;
    }

    template<ParallelPolicy P>
    [[nodiscard]] Task<> parallel_for(auto func, size_t num_loop) {
        assert(num_loop > 0);
        if constexpr (P == Concurrent) {
            for (size_t i = 0; i < num_loop; ++i)
                func(i);
            co_return;
        }
        else {
            Array<Task<Thread>> tasks(num_loop);
            for (size_t i = 0; i < num_loop; ++i) {
                tasks[i] = [](auto func, size_t i) noexcept -> Task<Thread> {
                    func(i);
                    co_return;
                }(func, i);
            }

            co_await std::suspend_always{};
            for (auto& task : tasks)
                task.wait();
        }
    }

    template<ParallelPolicy P>
    [[nodiscard]] Task<> parallel_for(auto func, size_t num_loop, int part) {
        assert(part > 0 && "[Error]: part must be a positive int");
        if constexpr (P == Concurrent) {
            for (size_t i = 0; i < num_loop; ++i)
                func(i);
            co_return;
        }
        else {
            using Range = std::pair<unsigned int, unsigned int>;
            static auto splitJob = [](size_t num_loop, int part, int i) {
                assert(0 <= i && i < part);
                const size_t maxLoopPerCore = (num_loop + part - 1) / part;
                const size_t from = i * maxLoopPerCore; 
                const size_t to = std::min(from + maxLoopPerCore, num_loop);
                return std::make_pair(from, to);
            };

            Array<Task<Thread>> tasks(part);
            for (int i = 0; i < part; ++i) {
                tasks[i] = [](auto func, Range range) noexcept -> Task<Thread> {
                    for (size_t loop = range.first; loop < range.second; ++loop)
                        func(loop);
                    co_return;
                }(func, splitJob(num_loop, part, i));
            }

            co_await std::suspend_always{};
            for (auto& task : tasks)
                task.wait();
        }
    }
}
