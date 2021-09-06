/*
 * Copyright 2021 WeiBo He.
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

namespace Physica::Core::Parallel {
    ThreadPool* ThreadPool::instance = nullptr;
    thread_local ThreadPool::ThreadInfo* ThreadPool::info = nullptr;

    ThreadPool::ThreadPool(unsigned int threadCount) : thread_data(threadCount), exit(false) {
        for (unsigned int i = 0; i < threadCount; ++i) {
            thread_data[i].thread.reset(new std::thread([this, i]() { workerMainLoop(i); } ));
        }
    }

    ThreadPool::~ThreadPool() {
        exit = true;
        for (auto& data : thread_data) {
            auto& thread = data.thread;
            if (thread->joinable())
                thread->join();
        }
    }

    void ThreadPool::workerMainLoop(unsigned int thread_id) {
        auto& threadInfo = getThreadInfo();
        threadInfo.pool = this;
        threadInfo.id = thread_id;
        auto& data = thread_data[thread_id];
        auto& queue = data.queue;
        std::unique_lock locker(data.queueMutex, std::defer_lock);
        while (true) {
            if (!queue.empty()) {
                locker.lock();
                std::unique_ptr<Task> task(std::move(queue.front()));
                queue.pop();
                locker.unlock();
                task->execute();
            }
            else if (exit)
                return;
            else
                std::this_thread::yield();
        }
    }

    ThreadPool::ThreadInfo& ThreadPool::getThreadInfo() {
        if (info == nullptr) {
            info = new ThreadInfo();
            info->pool = nullptr;
            info->randState = std::hash<std::thread::id>()(std::this_thread::get_id());
            info->id = 0;
        }
        return *info;
    }

    void ThreadPool::initThreadPool(unsigned int threadCount) {
        if (instance == nullptr) {
            instance = new ThreadPool(threadCount);
        }
    }
}
