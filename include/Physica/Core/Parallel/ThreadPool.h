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
#include <sys/sysinfo.h>

namespace Physica::Core::Parallel {
    /**
     * Reference:
     * Eigen https://eigen.tuxfamily.org/
     */
    class ThreadPool {
    public:
        struct ThreadData {
            std::unique_ptr<std::thread> thread;
            std::queue<std::function<void()>> queue;
            std::mutex queueMutex;

            ThreadData() : thread(), queue(), queueMutex() {}
        };

        struct ThreadInfo {
            ThreadPool* pool;
            uint64_t randState;
            size_t id;
        };
    private:
        thread_local static ThreadInfo* info;

        std::vector<ThreadData> thread_data;
        bool exit;
    public:
        ThreadPool(unsigned int threadCount);
        ThreadPool(const ThreadPool&) = delete;
        ThreadPool(ThreadPool&&) noexcept = delete;
        ~ThreadPool();
        /* Operators */
        ThreadPool& operator=(const ThreadPool&) = delete;
        ThreadPool& operator=(ThreadPool&&) noexcept = delete;
        /* Operations */
        void schedule(std::function<void()> func);
        void shouldExit() { exit = true; }
        /* Getters */
        [[nodiscard]] unsigned int getThreadCount() const noexcept { return thread_data.size(); }
        /* Static Members */
        [[nodiscard]] static ThreadInfo& getThreadInfo();
    private:
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
}