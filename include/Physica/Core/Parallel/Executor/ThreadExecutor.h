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

#include <type_traits>
#include <future>
#include <cassert>
#include "SeqExecutor.h"
#include "Physica/Core/Parallel/Task.h"

namespace Physica {
    class PHYSICA_API ThreadExecutor {
    public:
        using FutureType = std::future<void>;
        using Range = SeqExecutor::Range;
    public:
        /* Operations */
        template<class Functor>
        [[nodiscard("[Warn]: Discarding async return value or exception")]]
        static Task<Thread> schedule(Functor func) noexcept;
        template<class Functor>
        [[nodiscard]] static Task<> parallel_for(Functor func, size_t loopCount);
        template<class Functor>
        [[nodiscard]] static Task<> parallel_for(Functor func, size_t loopCount, int core);

        static void wait() {}
        /* Getters */
        [[nodiscard]] static int getThreadID() { return ThreadPool::getInstance().getThreadID(); }
        [[nodiscard]] static int getNumThread() { return ThreadPool::getInstance().getNumThreads(); }
        /* Static members */
        [[nodiscard]] inline static Range splitJob(size_t loopCount, int core, int part);
    };

    template<class Functor>
    Task<Thread> ThreadExecutor::schedule(Functor func) noexcept {
        func();
        co_return;
    }

    template<class Functor>
    Task<> ThreadExecutor::parallel_for(Functor func, size_t loopCount) {
        using ResultType = std::invoke_result<Functor, size_t>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        assert(loopCount > 0);
        Array<Task<Thread>> tasks(loopCount);
        for (size_t i = 0; i < loopCount; ++i) {
            tasks[i] = [](auto func, size_t i) noexcept -> Task<Thread> {
                func(i);
                co_return;
            }(func, i);
        }

        for (auto& task : tasks) {
            while (!task.done())
                co_await std::suspend_always{};
            task.get();
        }
    }

    template<class Functor>
    Task<> ThreadExecutor::parallel_for(Functor func, size_t loopCount, int core) {
        using ResultType = std::invoke_result<Functor, size_t>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        assert(core > 0 && "[Error]: core must be a positive int");
        Array<Task<Thread>> tasks(core);
        for (int i = 0; i < core; ++i) {
            tasks[i] = [](auto func, Range range) noexcept -> Task<Thread> {
                for (size_t loop = range.first; loop < range.second; ++loop)
                    func(loop);
                co_return;
            }(func, splitJob(loopCount, core, i));
        }

        for (auto& task : tasks) {
            while (!task.done())
                co_await std::suspend_always{};
            task.get();
        }
    }

    inline auto ThreadExecutor::splitJob(size_t loopCount, int core, int part) -> Range {
        assert(0 <= part && "[Error]: part must be a positive int");
        assert(part < core && "[Error]: More partition than core number requested");
        const size_t maxLoopPerCore = (loopCount + core - 1) / core;
        const size_t from = part * maxLoopPerCore; 
        const size_t to = std::min(from + maxLoopPerCore, loopCount);
        return std::make_pair(from, to);
    }
}

namespace Physica {
    template<>
    class Traits<ThreadExecutor> {
    public:
        constexpr static bool UseCPU = true;
        constexpr static bool UseCUDA = false;
    };
}
