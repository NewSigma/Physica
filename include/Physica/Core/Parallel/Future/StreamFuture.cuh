/*
 * Copyright 2023-2024 Weibo He.
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
        [[nodiscard]] static std::unique_ptr<StreamFuture> makeFuture(); //[Warn]: Discard a future may leak exceptions
    private:
        StreamFuture();
        static void taskDoneCallback(void* p_future);
    };
}
