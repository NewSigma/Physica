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
#include <Physica/Core/Parallel/Future/StreamFuture.cuh>
#include <Physica/Core/Parallel/CUDAContext.cuh>
#include <Physica/Core/Exception/CUDA/CUDA.cuh>

namespace Physica::Core {
    StreamFuture::StreamFuture() : isDone(false) {
        check(cudaLaunchHostFunc(CUDAContext::getInstance(), &StreamFuture::taskDoneCallback, this));
    }

    StreamFuture::~StreamFuture() {
        if (!isDone) {
            std::unique_lock locker(mutex);
            cond.wait(locker, [this]() { return isDone; });
        }
    }

    void StreamFuture::wait() {
        std::unique_lock locker(mutex);
        cond.wait(locker, [this]() { return isDone; });
        check(cudaGetLastError());
    }

    std::unique_ptr<StreamFuture> StreamFuture::makeFuture() {
        return std::unique_ptr<StreamFuture>(new StreamFuture());
    }

    void StreamFuture::taskDoneCallback(void* p_future) {
        auto& future = *reinterpret_cast<StreamFuture*>(p_future);
        std::unique_lock locker(future.mutex);
        future.isDone = true;
        locker.unlock();
        future.cond.notify_all();
    }
}
