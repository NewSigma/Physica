/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Parallel/Future/StreamFuture.cuh"
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"

namespace Physica {
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
