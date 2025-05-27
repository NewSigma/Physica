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
#include <utility>

namespace Physica {
    class SeqExecutor {
        class DummyFuture {
        public:
            constexpr static void wait() {}
            constexpr static void wait_async() {}
            constexpr static void get() {}
        };
    public:
        using Range = std::pair<unsigned int, unsigned int>;
    public:
        /* Operations */
        template<class Functor, class... Args>
        static DummyFuture schedule(Functor func, Args&&... args);
        template<class Functor>
        static DummyFuture parallel_for(Functor func, unsigned int loopCount);
        template<class Functor>
        static DummyFuture parallel_for(Functor func, unsigned int loopCount, unsigned int core);

        static void wait() {}
        /* Getters */
        [[nodiscard]] constexpr static unsigned int getNumThread() { return 1; }
        /* Static members */
        [[nodiscard]] inline static Range splitJob(unsigned int loopCount, unsigned int core, unsigned int part);
    };

    template<class Functor, class... Args>
    auto SeqExecutor::schedule(Functor func, Args&&... args) -> DummyFuture {
        func(std::forward<Args>(args)...);
        return {};
    }

    template<class Functor>
    auto SeqExecutor::parallel_for(Functor func, unsigned int loopCount) -> DummyFuture {
        using ResultType = std::invoke_result<Functor, unsigned int>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        assert(loopCount > 0);
        for (unsigned int i = 0; i < loopCount; ++i)
            func(i);
        return {};
    }

    template<class Functor>
    auto SeqExecutor::parallel_for(
            Functor func, unsigned int loopCount, unsigned int) -> DummyFuture {
        using ResultType = std::invoke_result<Functor, unsigned int>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        for (unsigned int i = 0; i < loopCount; ++i)
            func(i);
        return {};
    }

    inline auto SeqExecutor::splitJob(unsigned int loopCount, [[maybe_unused]] unsigned int core, [[maybe_unused]] unsigned int part) -> Range {
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
