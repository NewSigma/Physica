/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Parallel/ThreadPool.h"
#include <memory>
#ifdef PHYSICA_MKL
    #include <mkl_vml.h>
#endif

using namespace Physica;

namespace {


    void setThreadEnv() noexcept {
    #ifdef PHYSICA_MKL
        vmlSetMode(VML_HA | VML_FTZDAZ_CURRENT | VML_ERRMODE_NOERR);
    #endif
    }

    class GlobalEnv {
    public:
        GlobalEnv() noexcept {
            setThreadEnv();
        }
    };

    [[maybe_unused]] const GlobalEnv init{};

    int getNumProcesser() noexcept {
    #ifdef __linux__
        return get_nprocs();
    #else
        SYSTEM_INFO sinfo;
        GetSystemInfo(&sinfo);
        return sinfo.dwNumberOfProcessors;
    #endif
    }

    int makeNumThread() noexcept {
        const auto numProcesser = getNumProcesser();
        int num = ThreadPool::numThreadRequired;
        if (num == 0 || num > numProcesser)
            num = numProcesser * 3 / 4;
        return num;
    }
    /**
     * Generate a sequence of random seed (using the PCG-XSH-RS scheme)
     *
     * Reference:
     * [1] Eigen; https://eigen.tuxfamily.org  
     * [2] https://www.pcg-random.org/
     */
    uint32_t toNextState(uint64_t& state) noexcept {
        const uint64_t current = state;
        state = current * 6364136223846793005ULL + 0xda3e39cb94b95bdbULL;
        return static_cast<uint32_t>((current ^ (current >> 22U)) >> (22U + (current >> 61U)));
    }

    struct ThreadInfo {
        int id;
        uint64_t randState;
    };

    ThreadInfo& getThreadInfo() noexcept {
        thread_local static std::unique_ptr<ThreadInfo> info = nullptr;
        if (info == nullptr) {
            info = std::make_unique<ThreadInfo>();
            info->id = ThreadPool::MainThreadID;
            info->randState = std::hash<std::thread::id>()(std::this_thread::get_id());
        }
        return *info;
    }
}

int ThreadPool::numThreadRequired = 0;

void ThreadPool::ThreadData::push(Handle handle) noexcept {
    std::unique_lock locker(mutex);
    queue.push(handle);
}

auto ThreadPool::ThreadData::pop() noexcept -> Handle {
    Handle handle = nullptr;
    std::unique_lock locker(mutex);
    if (!queue.empty()) {
        handle = queue.front();
        queue.pop();
    }
    return handle;
}

ThreadPool::ThreadPool(int numThreads) : thread_data(numThreads), exit(false) {
    assert(numThreads > 0 && "[Error]: numThreads must be positive");
    for (int i = 0; i < numThreads; ++i) {
        thread_data[i].thread = std::thread([i]() noexcept {
            setThreadEnv();
            getInstance().workerMainLoop(i);
        });
    }
}

ThreadPool::~ThreadPool() {
    waitExit();
}

void ThreadPool::schedule(Handle handle) noexcept {
    assert(handle != nullptr);
    assert(!handle.done());
    int schedule_to{};
    if (isMainThread()) {
        static int MainThreadCounter = 0;
        schedule_to = MainThreadCounter;
        MainThreadCounter = (MainThreadCounter + 1) % getNumThreads();
    }
    else
        schedule_to = getThreadInfo().id;

    thread_data[schedule_to].push(handle);
    cond.notify_one();
}

auto ThreadPool::steal() noexcept -> Handle {
    const auto random = toNextState(getThreadInfo().randState);
    const int numThreads = getNumThreads();
    for (int i = 0; i < numThreads; ++i) {
        Handle handle = thread_data[(random + i) % numThreads].pop();
        if (handle)
            return handle;
    }
    return nullptr;
}

void ThreadPool::waitExit() {
    shouldExit();
    for (auto& data : thread_data) {
        auto& thread = data.thread;
        if (thread.joinable())
            thread.join();
    }
}

void ThreadPool::restart() {
    waitExit();
    exit = false;
    const int numThread = makeNumThread();
    thread_data = Array<ThreadData>(numThread);
    for (int i = 0; i < numThread; ++i) {
        thread_data[i].thread = std::thread([i]() noexcept {
            getInstance().workerMainLoop(i);
        });
    }
}

void ThreadPool::shouldExit() noexcept {
    poolMutex.lock();
    exit = true;
    poolMutex.unlock();
    cond.notify_all();
}

auto ThreadPool::getInstance() noexcept -> This& {
    static ThreadPool pool(makeNumThread());
    return pool;
}

int ThreadPool::getThreadID() noexcept {
    return getThreadInfo().id;
}

void ThreadPool::workerMainLoop(int thread_id) noexcept {
    getThreadInfo().id = thread_id;
    auto& data = thread_data[thread_id];
    while (true) {
        Handle handle = data.pop();
        if (!handle)
            handle = steal();

        if (handle) {
            assert(!handle.done() && "Data race");
            cond.notify_one();
            handle.resume();
        }
        else {
            std::unique_lock poolLocker(poolMutex);
            if (exit)
                return;
            cond.wait(poolLocker);
        }
    }
}
