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
#pragma once

#include <memory>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <vector>
#include <future>
#include <sys/sysinfo.h>
#include "PackagedTask.h"

namespace Physica::Core::Parallel {
    /**
     * Reference:
     * Eigen https://eigen.tuxfamily.org/
     */
    class ThreadPool {
    public:
        struct ThreadData {
            std::unique_ptr<std::thread> thread;
            std::queue<std::unique_ptr<Task>> queue;
            std::mutex queueMutex;

            ThreadData() : thread(), queue(), queueMutex() {}
        };

        struct ThreadInfo {
            ThreadPool* pool;
            uint64_t randState;
            size_t id;
        };
    private:
        static ThreadPool* instance;
        thread_local static ThreadInfo* info;

        std::vector<ThreadData> thread_data;
        bool exit;
    public:
        ThreadPool(const ThreadPool&) = delete;
        ThreadPool(ThreadPool&&) noexcept = delete;
        ~ThreadPool();
        /* Operators */
        ThreadPool& operator=(const ThreadPool&) = delete;
        ThreadPool& operator=(ThreadPool&&) noexcept = delete;
        /* Operations */
        template<class Function, class... Args>
        std::future<typename std::invoke_result<Function>::type> schedule(Function func, Args... args);
        void shouldExit() { exit = true; }
        /* Getters */
        [[nodiscard]] unsigned int getThreadCount() const noexcept { return thread_data.size(); }
        /* Static Members */
        [[nodiscard]] static ThreadInfo& getThreadInfo();
        static void initThreadPool(unsigned int threadCount);
        static void deInitThreadPool() { delete instance; instance = nullptr; }
        [[nodiscard]] static ThreadPool& getInstance() { return *instance; }
    private:
        ThreadPool(unsigned int threadCount);
        void workerMainLoop(unsigned int thread_id);
        /* Static Members */
        [[nodiscard]] static inline unsigned int defaultThreadNum() noexcept { return get_nprocs() * 3 / 4; }
        [[nodiscard]] static inline unsigned int threadRand(uint64_t& state);
    };

    inline unsigned int ThreadPool::threadRand(uint64_t& state) {
        uint64_t current = state;
        state = current * 6364136223846793005ULL + 0xda3e39cb94b95bdbULL;
        // Generate the random output (using the PCG-XSH-RS scheme)
        return static_cast<unsigned int>((current ^ (current >> 22U)) >> (22U + (current >> 61U)));
    }

    template<class Function, class ... Args>
    std::future<typename std::invoke_result<Function>::type> ThreadPool::schedule(Function func, Args... args) {
        using ResultType = typename std::invoke_result<Function>::type;
        std::packaged_task<ResultType(Args...)> task(std::bind(func, args...));
        auto result = task.get_future();
        unsigned int random_id = threadRand(getThreadInfo().randState) % thread_data.size();
        auto& data = thread_data[random_id];
        std::lock_guard lock(data.queueMutex);
        data.queue.emplace(new PackagedTask(std::move(task)));
        return result;
    }
}
