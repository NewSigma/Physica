/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

using namespace Physica;

void ThreadExecutor::auto_wait(FutureType& future) {
    const auto nano = std::chrono::nanoseconds(1);
    while (future.wait_for(nano) != std::future_status::ready) {
        std::unique_ptr<Task> task = ThreadPool::getInstance().steal();
        if (task != nullptr)
            task->execute();
        else
            std::this_thread::yield();
    }
    future.get();
}
