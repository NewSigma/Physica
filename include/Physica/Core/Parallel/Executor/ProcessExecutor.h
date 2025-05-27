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

#include <cassert>
#include "Physica/Core/Utils/Container/Array.h"
#include "Physica/Core/Parallel/Future/ProcessFuture.h"
#include "Physica/Core/Parallel/SubProcess.h"
#include "Physica/Core/Parallel/Parallel.h"

namespace Physica {
    class ProcessExecutor {
    public:
        using FutureType = ProcessFuture;
    public:
        PHYSICA_API static int nice_incr;
    public:
        /* Operations */
        template<class Functor, class... Args>
        static FutureType schedule(Functor func, Args&&... args);
        template<class Functor>
        static Task<Concurrent> parallel_for(Functor func, unsigned int loopCount, unsigned int core);
    };

    template<class Functor, class... Args>
    auto ProcessExecutor::schedule(Functor func, Args&&... args) -> FutureType {
        using ResultType = std::invoke_result<Functor, Args&&...>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: ProcessExecutor does not support functors with return value");

        SubProcess process([=]() { func(std::forward<Args>(args)...); }, nice_incr);
        return process.execute();
    }

    template<class Functor>
    Task<Concurrent> ProcessExecutor::parallel_for(Functor func, unsigned int loopCount, unsigned int core) {
        using ResultType = std::invoke_result<Functor, unsigned int>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        assert(loopCount >= core);
        assert(core > 0);

        const unsigned int maxLoopPerCore = (loopCount + core - 1) / core;
        unsigned int from = 0; 
        unsigned int to = maxLoopPerCore;
        Array<FutureType> futures(core);
        for (unsigned int _ = 0; _ < core; ++_) {
            SubProcess process([=]() {
                                   for (unsigned int i = from; i < to; ++i)
                                       func(i);
                               }, nice_incr);
            futures.append(process.execute());
            from += maxLoopPerCore;
            const unsigned int next_to = to + maxLoopPerCore;
            to = next_to > loopCount ? loopCount : next_to;
        }
        co_await std::suspend_always{};

        for (auto& future : futures)
            std::ignore = future.wait();
    }
}
