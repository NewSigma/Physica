/*
 * Copyright 2022-2024 Weibo He.
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
#include <Physica/Core/Parallel/Future/FutureGroup.h>
#include <Physica/Core/Parallel/ThreadPool.h>
#include "SequentialExecutor.h"

namespace Physica::Core {
    class ThreadExecutor {
    public:
        using FutureType = std::future<void>;
        using Range = typename SequentialExecutor::Range;
    public:
        /* Operations */
        template<class Functor, class... Args>
        [[nodiscard]] static FutureType schedule(Functor func, Args&&... args) noexcept;
        template<class Functor>
        [[nodiscard]] static FutureGroup<FutureType> parallel_for(Functor func, unsigned int loopCount);
        template<class Functor>
        [[nodiscard]] static FutureGroup<FutureType> parallel_for(Functor func, unsigned int loopCount, unsigned int core);
        PHYSICA_API static void auto_wait(FutureType& future);
        inline static void auto_wait(FutureGroup<FutureType>& group);
        static void wait() {}
        /* Getters */
        [[nodiscard]] static unsigned int getNumThread() { return ThreadPool::getInstance().getNumThreads(); }
        /* Static members */
        [[nodiscard]] inline static Range splitJob(unsigned int loopCount, unsigned int core, unsigned int part);
    };
    /**
     * Return value are declared as [[nodiscard]] to warn user pay attention to
     * 1. asynchronous return value
     * 2. exception stored in std::future, which is usually ignored
     */
    template<class Functor, class... Args>
    typename ThreadExecutor::FutureType ThreadExecutor::schedule(Functor func, Args&&... args) noexcept {
        return ThreadPool::getInstance().schedule(std::move(func), std::forward<Args>(args)...);
    }

    template<class Functor>
    FutureGroup<typename ThreadExecutor::FutureType> ThreadExecutor::parallel_for(
            Functor func, unsigned int loopCount) {
        using ResultType = typename std::invoke_result<Functor, unsigned int>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        assert(loopCount > 0);
        FutureGroup<FutureType> result(loopCount);
        for (unsigned int i = 0; i < loopCount; ++i) {
            constexpr bool isNothrow = std::is_nothrow_invocable<Functor, unsigned int>::value;
            result.append(ThreadPool::getInstance().schedule([func, i]() noexcept(isNothrow) -> void { func(i); }));
        }
        return result;
    }

    template<class Functor>
    FutureGroup<typename ThreadExecutor::FutureType> ThreadExecutor::parallel_for(
            Functor func, unsigned int loopCount, unsigned int core) {
        using ResultType = typename std::invoke_result<Functor, unsigned int>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        assert(core > 0);
        FutureGroup<FutureType> result(core);
        for (unsigned int i = 0; i < core; ++i) {
            constexpr bool isNothrow = std::is_nothrow_invocable<Functor, unsigned int>::value;
            const auto range = splitJob(loopCount, core, i);
            result.append(ThreadPool::getInstance().schedule([range, func]() noexcept(isNothrow) -> void {
                for (unsigned int loop = range.first; loop < range.second; ++loop)
                    func(loop);
            }));
        }
        return result;
    }

    inline void ThreadExecutor::auto_wait(FutureGroup<FutureType>& group) {
        for (auto& future : group.getFutures())
            auto_wait(future);
    }

    inline typename ThreadExecutor::Range ThreadExecutor::splitJob(unsigned int loopCount, unsigned int core, unsigned int part) {
        assert(core != 0 && "[Error]: Zero core is not allowed");
        assert(part < core && "[Error]: More partition than core number requested");
        const unsigned int maxLoopPerCore = (loopCount + core - 1) / core;
        const unsigned int from = part * maxLoopPerCore; 
        const unsigned int to = std::min(from + maxLoopPerCore, loopCount);
        return std::make_pair(from, to);
    }
}

namespace Physica {
    template<>
    class Traits<Core::ThreadExecutor> {
    public:
        constexpr static bool UseCPU = true;
        constexpr static bool UseCUDA = false;
    };
}
