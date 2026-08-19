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

#include "Physica/Core/Parallel/Parallel.h"
#include "Physica/Core/Parallel/Task/Task.h"
#include "Physica/Core/Parallel/Task/When.h"
#include "Physica/Core/Parallel/ThreadPool.h"

namespace Physica {
    template<ExecutePolicy P>
    [[nodiscard]] Task schedule(std::invocable<> auto fn) requires(P == Thread) {
        co_await ThreadPool::getInstance();
        fn();
    }

    template<ExecutePolicy P>
    [[nodiscard]] Task parallel_for(std::invocable<size_t> auto fn, size_t num_loop) noexcept requires(P == Thread) {
        auto tasks = Array<Task>::generate([fn](size_t i) {
            return [](auto fn, size_t i) static noexcept -> Task {
                co_await toAwaiter(ThreadPool::getInstance()).implicit();
                fn(i);
            }(fn, i);
        }, num_loop);
        ThreadPool::getInstance().notify_one();
        co_await when_all(tasks);
    }

    template<ExecutePolicy P>
    [[nodiscard]] Task parallel_for(std::invocable<size_t> auto fn, size_t num_loop, size_t part) noexcept requires(P == Thread) {
        if (part == 0 || part > num_loop)
            part = std::min<size_t>(ThreadPool::getInstance().getNumThreads(), num_loop);
        assume(part > 0);

        auto tasks = Array<Task>::generate([fn, num_loop, part](int i) noexcept {
            using Range = std::pair<size_t, size_t>;
            const size_t maxLoopPerCore = (num_loop + part - 1) / part;
            const size_t from = i * maxLoopPerCore;
            const size_t to = std::min(from + maxLoopPerCore, num_loop);
            return [](auto fn, Range range) static noexcept -> Task {
                co_await toAwaiter(ThreadPool::getInstance()).implicit();
                for (size_t loop = range.first; loop < range.second; ++loop)
                    fn(loop);
            }(fn, std::make_pair(from, to));
        }, part);
        ThreadPool::getInstance().notify_one();
        co_await when_all(tasks);
    }
}
