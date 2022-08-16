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
#include <cassert>
#include "Physica/Core/Parallel/FutureGroup.h"

namespace Physica::Core::Parallel {
    class SequentialExecutor {
        using FutureType = DummyFuture;
    public:
        /* Operations */
        template<class Functor, class... Args>
        static FutureType schedule(Functor func, Args... args);
        template<class Functor>
        static FutureGroup<FutureType> parallel_for(Functor func, unsigned int loopCount, [[maybe_unused]] unsigned int core);
        /* Getters */
        [[nodiscard]] constexpr static unsigned int getNumThread() { return 1; }
    };

    template<class Functor, class... Args>
    typename SequentialExecutor::FutureType SequentialExecutor::schedule(Functor func, Args... args) {
        func(std::forward<Args>(args)...);
        return FutureType{};
    }

    template<class Functor>
    FutureGroup<typename SequentialExecutor::FutureType> SequentialExecutor::parallel_for(
            Functor func, unsigned int loopCount, [[maybe_unused]] unsigned int core) {
        using ResultType = typename std::invoke_result<Functor, unsigned int>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: Invalid functor");
        assert(loopCount >= core);
        assert(core > 0);
        for (unsigned int i = 0; i < loopCount; ++i)
            func(i);
        return FutureGroup<FutureType>{};
    }
}
