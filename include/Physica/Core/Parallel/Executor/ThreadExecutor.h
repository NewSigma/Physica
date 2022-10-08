/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Core/Parallel/Future/FutureGroup.h"
#include "Physica/Core/Parallel/ThreadPool.h"

namespace Physica::Core::Parallel {
    class ThreadExecutor {
    public:
        using FutureType = std::future<void>;
    public:
        /* Operations */
        template<class Functor, class... Args>
        static FutureType schedule(Functor func, Args... args);
        template<class Functor>
        static FutureGroup<FutureType> parallel_for(Functor func, unsigned int loopCount, unsigned int core);
        static void wait(FutureType& future);
        /* Getters */
        [[nodiscard]] static unsigned int getNumThread() { return ThreadPool::getInstance().getThreadCount(); }
    };

    template<class Functor, class... Args>
    typename ThreadExecutor::FutureType ThreadExecutor::schedule(Functor func, Args... args) {
        return ThreadPool::getInstance().schedule(std::move(func), std::forward<Args>(args)...);
    }

    template<class Functor>
    FutureGroup<typename ThreadExecutor::FutureType> ThreadExecutor::parallel_for(
            Functor func, unsigned int loopCount, unsigned int core) {
        using ResultType = typename std::invoke_result<Functor, unsigned int>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        assert(loopCount >= core);
        assert(core > 0);
        const unsigned int maxLoopPerCore = (loopCount + core - 1) / core;
        unsigned int from = 0; 
        unsigned int to = maxLoopPerCore;
        FutureGroup<FutureType> result(core);
        for (unsigned int _ = 0; _ < core; ++_) {
            result.append(ThreadPool::getInstance().schedule(
                [=]() -> void {
                    for (unsigned int i = from; i < to; ++i)
                        func(i);
                }));
            from += maxLoopPerCore;
            const unsigned int next_to = to + maxLoopPerCore;
            to = next_to > loopCount ? loopCount : next_to;
        }
        return result;
    }
}
