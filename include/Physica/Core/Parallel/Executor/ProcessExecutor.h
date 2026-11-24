/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
