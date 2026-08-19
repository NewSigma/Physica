/*
 * Copyright 2021-2026 Weibo He.
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
#include <xmmintrin.h>
#ifdef PHYSICA_MKL
    #include <mkl_vml.h>
#endif
#ifdef PHYSICA_HDF5
    #include <H5Epublic.h>
#endif

using namespace Physica;

namespace {
    void setThreadEnv() noexcept {
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    #ifdef PHYSICA_MKL
        vmlSetMode(VML_HA | VML_FTZDAZ_CURRENT | VML_ERRMODE_NOERR);
    #endif
    }

    class GlobalEnv {
    public:
        GlobalEnv() noexcept {
            setThreadEnv();
        #ifdef PHYSICA_HDF5
            H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
        #endif
        }
    };

    [[maybe_unused]] const GlobalEnv init{};

    int getNumProcessor() noexcept {
    #ifdef __linux__
        return get_nprocs();
    #else
        SYSTEM_INFO sinfo;
        GetSystemInfo(&sinfo);
        return sinfo.dwNumberOfProcessors;
    #endif
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

    auto& getThreadInfo() noexcept {
        thread_local static struct ThreadInfo {
            int id = ThreadPool::MainThreadID;
            uint64_t randState = std::hash<std::thread::id>()(std::this_thread::get_id());
        } info{};
        return info;
    }
}

class ThreadPool::ThreadQueue {
public:
    std::thread thread;
private:
    std::queue<Handle> queue;
    std::mutex mutex;
public:
    ThreadQueue() = default;
    /* Operations */
    void push(Handle handle) noexcept;
    [[nodiscard]] Handle pop() noexcept;
};

void ThreadPool::ThreadQueue::push(Handle handle) noexcept {
    std::lock_guard locker(mutex);
    queue.push(handle);
}

auto ThreadPool::ThreadQueue::pop() noexcept -> Handle {
    Handle handle = nullptr;
    std::lock_guard locker(mutex);
    if (!queue.empty()) {
        handle = queue.front();
        queue.pop();
    }
    return handle;
}

void ThreadPool::Scheduler::await_suspend(Handle handle) noexcept {
    assert(handle != nullptr);
    assert(!handle.done());
    int schedule_to{};
    auto& pool = getInstance();
    if (isMainThread()) {
        static int MainThreadCounter = 0;
        schedule_to = MainThreadCounter;
        MainThreadCounter = (MainThreadCounter + 1) % pool.getNumThreads();
    }
    else
        schedule_to = getThreadInfo().id;

    pool.queues[schedule_to].push(handle);
    pool.cond.notify_one();
}

void ThreadPool::Scheduler::on_wait(Handle) noexcept {
    if (auto handle = getInstance().steal())
        handle.resume();
    else
        std::this_thread::yield();
}

int ThreadPool::numThreadRequired = 0;

ThreadPool::ThreadPool(int numThreads) : queues(numThreads) {
    assert(numThreads > 0 && "[Error]: numThreads must be positive");
    for (int i = 0; i < numThreads; ++i)
        queues[i].thread = std::thread(&ThreadPool::workerMainLoop, this, i);
}

ThreadPool::~ThreadPool() {
    waitExit();
}

auto ThreadPool::operator co_await() noexcept -> Scheduler {
    return Scheduler{};
}

void ThreadPool::notify_one() {
    cond.notify_one();
}

void ThreadPool::notify_all() {
    cond.notify_all();
}

auto ThreadPool::steal() noexcept -> Handle {
    const auto random = toNextState(getThreadInfo().randState);
    const int numThreads = getNumThreads();
    for (int i = 0; i < numThreads; ++i) {
        Handle handle = queues[(random + i) % numThreads].pop();
        if (handle)
            return handle;
    }
    return nullptr;
}

void ThreadPool::shouldExit() noexcept {
    {
        auto locker = std::lock_guard(poolMutex);
        exit = true;
    }
    cond.notify_all();
}

void ThreadPool::waitExit() {
    shouldExit();
    for (auto& queue : queues) {
        auto& thread = queue.thread;
        if (thread.joinable())
            thread.join();
    }
}

void ThreadPool::restart() {
    waitExit();
    exit = false;
    const int numThread = getNumThreads();
    queues = Array<ThreadQueue>(numThread);
    for (int i = 0; i < numThread; ++i)
        queues[i].thread = std::thread(&ThreadPool::workerMainLoop, this, i);
}

auto ThreadPool::getInstance() noexcept -> This& {
    static ThreadPool pool([]() {
        int numProcessor = getNumProcessor();
        int num = ThreadPool::numThreadRequired;
        if (num == 0 || num > numProcessor)
            num = std::max(1, numProcessor * 3 / 4);
        return num;
    }());
    return pool;
}

int ThreadPool::getThreadID() noexcept {
    return getThreadInfo().id;
}

void ThreadPool::spin() noexcept { __builtin_ia32_pause(); }

void ThreadPool::workerMainLoop(int thread_id) noexcept {
    setThreadEnv();
    getThreadInfo().id = thread_id;
    auto& queue = queues[thread_id];
    while (true) {
        Handle handle = queue.pop();
        if (!handle)
            handle = steal();

        if (handle) {
            assert(!handle.done() && "Data race");
            cond.notify_one();
            handle.resume();
        }
        else {
            std::unique_lock locker(poolMutex);
            if (exit)
                return;
            cond.wait(locker);
        }
    }
}
