/*
 * Copyright 2021-2024 Weibo He.
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
#include <Physica/Core/Parallel/ThreadPool.h>

namespace Physica::Core {
    unsigned int ThreadPool::numThreadRequired = 0;
    thread_local std::unique_ptr<ThreadPool::ThreadInfo> ThreadPool::info = nullptr;

    ThreadPool::ThreadPool(unsigned int numThreads) : thread_data(numThreads), exit(false) {
        for (unsigned int i = 0; i < numThreads; ++i) {
            thread_data[i].thread.reset(new std::thread([this, i]() noexcept { workerMainLoop(i); } ));
        }
    }

    ThreadPool::~ThreadPool() {
        waitExit();
    }

    std::unique_ptr<Task> ThreadPool::steal() {
        const unsigned int numThreads = getNumThreads();
        const unsigned int random = threadRand(getThreadInfo().randState);
        for (size_t i = 0; i < numThreads; ++i) {
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
        const unsigned int numThread = makeNumThread();
        thread_data = Utils::Array<ThreadData>(numThread);
        for (unsigned int i = 0; i < numThread; ++i) {
            thread_data[i].thread.reset(new std::thread([this, i]() { workerMainLoop(i); } ));
        }
    }

    void ThreadPool::shouldExit() noexcept {
        exit = true;
        cond.notify_all();
    }

    ThreadPool::ThreadInfo& ThreadPool::getThreadInfo() noexcept {
        if (info == nullptr) {
            info.reset(new ThreadInfo());
            info->id = MainThreadID;
            info->numScheduled = 0;
            info->randState = std::hash<std::thread::id>()(std::this_thread::get_id());
        }
        return *info;
    }

    ThreadPool& ThreadPool::getInstance() noexcept {
        static ThreadPool pool(makeNumThread());
        return pool;
    }

    void ThreadPool::workerMainLoop(unsigned int thread_id) noexcept {
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
            else if (!exit) {
                std::unique_lock poolLocker(poolMutex);
                cond.wait(poolLocker);
            }
            else
                return;
        }
    }
}
