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
#include "Physica/Core/Math/Random/RandomSeed.h"

using namespace Physica;

int ThreadPool::numThreadRequired = 0;

ThreadPool::ThreadPool(int numThreads)
        : thread_data(numThreads), exit(false) {
    assert(numThreads > 0 && "[Error]: numThreads must be positive");
    for (int i = 0; i < numThreads; ++i) {
        thread_data[i].thread.reset(new std::thread([this, i]() noexcept { workerMainLoop(i); }));
    }
}

ThreadPool::~ThreadPool() {
    waitExit();
}

std::unique_ptr<Task> ThreadPool::steal() {
    const auto random = RandomSeed::toNextSeed(getThreadInfo().randState);
    const int numThreads = getNumThreads();
    for (int i = 0; i < numThreads; ++i) {
        ThreadData& data = thread_data[(random + i) % numThreads];
        std::unique_lock locker(data.queueMutex);
        auto& queue = data.queue;
        if (!queue.empty()) {
            std::unique_ptr<Task> task(std::move(queue.front()));
            queue.pop();
            return task;
        }
    }
    return std::unique_ptr<Task>(nullptr);
}

void ThreadPool::waitExit() {
    shouldExit();
    for (auto& data : thread_data) {
        auto& thread = data.thread;
        if (thread->joinable())
            thread->join();
    }
}

void ThreadPool::restart() {
    waitExit();
    exit = false;
    const int numThread = makeNumThread();
    thread_data = Array<ThreadData>(numThread);
    for (int i = 0; i < numThread; ++i) {
        thread_data[i].thread.reset(new std::thread([this, i]() { workerMainLoop(i); }));
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

void ThreadPool::workerMainLoop(int thread_id) noexcept {
    auto& threadInfo = getThreadInfo();
    threadInfo.id = thread_id;
    auto& data = thread_data[thread_id];
    auto& queue = data.queue;
    std::unique_lock locker(data.queueMutex, std::defer_lock);
    while (true) {
        std::unique_ptr<Task> task = nullptr;
        locker.lock();
        if (!queue.empty()) {
            task = std::move(queue.front());
            queue.pop();
            locker.unlock();
        }
        else {
            locker.unlock();
            task = steal();
        }

        if (task != nullptr) {
            cond.notify_one();
            task->execute();
        }
        else {
            std::unique_lock poolLocker(poolMutex);
            if (exit)
                return;
            cond.wait(poolLocker);
        }
    }
}

auto ThreadPool::getThreadInfo() noexcept -> ThreadInfo& {
    thread_local static std::unique_ptr<ThreadInfo> info = nullptr;
    if (info == nullptr) {
        info.reset(new ThreadInfo());
        info->id = MainThreadID;
        info->randState = std::hash<std::thread::id>()(std::this_thread::get_id());
    }
    return *info;
}

int ThreadPool::getNumProcesser() noexcept {
#ifdef __linux__
    return get_nprocs();
#else
    SYSTEM_INFO sinfo;
    GetSystemInfo(&sinfo);
    return sinfo.dwNumberOfProcessors;
#endif
}

int ThreadPool::makeNumThread() noexcept {
    const auto numProcesser = getNumProcesser();
    if (numThreadRequired == 0 || numThreadRequired > numProcesser)
        numThreadRequired = numProcesser * 3 / 4;
    return numThreadRequired;
}
