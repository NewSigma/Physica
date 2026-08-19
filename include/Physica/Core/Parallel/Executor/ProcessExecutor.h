/*
 * Copyright 2022-2026 Weibo He.
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
#include "Physica/Core/Parallel/Future/ProcessFuture.h"
#include "Physica/Core/Parallel/SubProcess.h"

namespace Physica {
    class ProcessExecutor {
    public:
        using FutureType = ProcessFuture;
    public:
        PHYSICA_API static int nice_incr;
    public:
        /* Operations */
        template<class... Args>
        static FutureType schedule(std::invocable<Args...> auto fn, Args&&... args);
    };

    template<class... Args>
    auto ProcessExecutor::schedule(std::invocable<Args...> auto fn, Args&&... args) -> FutureType {
        using ResultType = std::invoke_result<decltype(fn), Args&&...>::type;
        static_assert(std::is_same<void, ResultType>::value, "[Error]: ProcessExecutor does not support functors with return value");

        SubProcess process([=]() { fn(std::forward<Args>(args)...); }, nice_incr);
        return process.execute();
    }
}
