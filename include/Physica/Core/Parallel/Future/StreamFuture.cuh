/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <mutex>
#include <condition_variable>
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API StreamFuture {
        bool isDone;
        std::mutex mutex;
        std::condition_variable cond;
    public:
        StreamFuture(const StreamFuture&) = delete;
        StreamFuture(StreamFuture&&) noexcept = delete;
        ~StreamFuture();
        /* Operators */
        StreamFuture& operator=(const StreamFuture&) = delete;
        StreamFuture& operator=(StreamFuture&&) noexcept = delete;
        /* Operations */
        void wait();
        /* Static members */
        [[nodiscard("[Warn]: Discard a future may leak exceptions")]]
        static std::unique_ptr<StreamFuture> makeFuture(); //
    private:
        StreamFuture();
        static void taskDoneCallback(void* p_future);
    };
}
