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
#include <cassert>
#include "Physica/Core/Parallel/Future/DummyFuture.h"

namespace Physica {
    class SeqExecutor {
    public:
        using FutureType = DummyFuture;
        using Range = std::pair<unsigned int, unsigned int>;
    public:
        /* Operations */
        template<class Functor, class... Args>
        inline static FutureType schedule(Functor func, Args&&... args);
        template<class Functor>
        inline static FutureGroup<FutureType> parallel_for(Functor func, unsigned int loopCount);
        template<class Functor>
        inline static FutureGroup<FutureType> parallel_for(Functor func, unsigned int loopCount, unsigned int core);
        static void auto_wait(FutureType&) {}
        static void auto_wait(FutureGroup<FutureType>&) {}
        static void wait() {}
        /* Getters */
        [[nodiscard]] constexpr static unsigned int getNumThread() { return 1; }
        /* Static members */
        [[nodiscard]] inline static Range splitJob(unsigned int loopCount, unsigned int core, unsigned int part);
    };

    template<class Functor, class... Args>
    inline SeqExecutor::FutureType SeqExecutor::schedule(Functor func, Args&&... args) {
        func(std::forward<Args>(args)...);
        return FutureType{};
    }

    template<class Functor>
    inline FutureGroup<typename SeqExecutor::FutureType> SeqExecutor::parallel_for(
            Functor func, unsigned int loopCount) {
        using ResultType = std::invoke_result<Functor, unsigned int>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        assert(loopCount > 0);
        for (unsigned int i = 0; i < loopCount; ++i)
            func(i);
        return FutureGroup<FutureType>{};
    }

    template<class Functor>
    inline FutureGroup<typename SeqExecutor::FutureType> SeqExecutor::parallel_for(
            Functor func, unsigned int loopCount, unsigned int) {
        using ResultType = std::invoke_result<Functor, unsigned int>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        for (unsigned int i = 0; i < loopCount; ++i)
            func(i);
        return FutureGroup<FutureType>{};
    }

    inline SeqExecutor::Range SeqExecutor::splitJob(
            unsigned int loopCount, [[maybe_unused]] unsigned int core, [[maybe_unused]] unsigned int part) {
        assert(part == 0 && core == 1 && "[Error]: SeqExecutor do not support other param");
        return std::make_pair(0, loopCount);
    }
}

namespace Physica {
    template<class T> class Traits;

    template<>
    class Traits<SeqExecutor> {
    public:
        constexpr static bool UseCPU = true;
        constexpr static bool UseCUDA = false;
    };
}
