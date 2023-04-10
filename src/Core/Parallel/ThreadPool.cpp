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
        cond.notify_all();
        for (auto& data : thread_data) {
            auto& thread = data.thread;
            if (thread->joinable())
                thread->join();
        }
    }

    std::unique_ptr<Task> ThreadPool::steal() {
        const unsigned int threadCount = getThreadCount();
        const unsigned int random = threadRand(getThreadInfo().randState);
        for (size_t i = 0; i < threadCount; ++i) {
            ThreadData& data = thread_data[(random + i) % threadCount];
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

    void ThreadPool::workerMainLoop(unsigned int thread_id) {
        bindToCore(thread_id);
        auto& threadInfo = getThreadInfo();
        threadInfo.id = thread_id;
        auto& data = thread_data[thread_id];
        auto& queue = data.queue;
        std::unique_lock locker(data.queueMutex, std::defer_lock);
        while (true) {
            locker.lock();
            if (!queue.empty()) {
                std::unique_ptr<Task> task(std::move(queue.front()));
                queue.pop();
                locker.unlock();
                task->execute();
            }
            else {
                locker.unlock();
                std::unique_ptr<Task> task = steal();
                if (task != nullptr)
                    task->execute();
                else if (exit)
                    return;
                else {
                    std::unique_lock poolLocker(poolMutex);
                    cond.wait(poolLocker);
                }
            }
        }
    }

    void ThreadPool::bindToCore(unsigned int thread_id) {
        cpu_set_t set;
        CPU_ZERO(&set);
        CPU_SET(thread_id, &set);
        const pthread_t thread = thread_data[thread_id].thread->native_handle();
        pthread_setaffinity_np(thread, sizeof(cpu_set_t), &set);
    }

    ThreadPool::ThreadInfo& ThreadPool::getThreadInfo() {
        if (info == nullptr) {
            info = new ThreadInfo();
            info->id = MainThreadID;
            info->numScheduled = 0;
            info->randState = std::hash<std::thread::id>()(std::this_thread::get_id());
        }
        return *info;
    }

    void ThreadPool::initThreadPool(unsigned int threadCount) {
        if (instance == nullptr) {
            instance = new ThreadPool(threadCount);
        }
    }
}
