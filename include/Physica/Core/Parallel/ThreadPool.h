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
#pragma once

#ifdef __linux__
    #include <sys/sysinfo.h>
#else
    #include <windows.h>
#endif
#include <memory>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <limits>
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    /**
     * Reference:
     * [1] Eigen; https://eigen.tuxfamily.org
     */
    class PHYSICA_API ThreadPool final {
        using This = ThreadPool;
        using Handle = std::coroutine_handle<>;
        struct ThreadData {
            std::unique_ptr<std::thread> thread;
            std::queue<Handle> queue;
            std::mutex queueMutex;

            ThreadData() : thread(), queue(), queueMutex() {}
        };

        struct ThreadInfo {
            int id;
            uint64_t randState;
        };
    public:
        static int numThreadRequired;
    private:
        // Larger than any thread ID, main thread is not maintained by thread pool
        constexpr static int MainThreadID = std::numeric_limits<decltype(numThreadRequired)>::max();

        Array<ThreadData> thread_data;
        std::mutex poolMutex;
        std::condition_variable cond;
        bool exit;
    public:
        ThreadPool(const This&) = delete;
        ThreadPool(This&&) noexcept = delete;
        ~ThreadPool();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void schedule(Handle handle) noexcept;
        [[nodiscard]] Handle steal();
        void waitExit();
        void restart();
        /* Getters */
        [[nodiscard]] int getNumThreads() const noexcept { return thread_data.getLength(); }
        /* Setters */
        void shouldExit() noexcept;
        /* Static Members */
        [[nodiscard]] static This& getInstance() noexcept;
        [[nodiscard]] static inline int getThreadID() noexcept;
        [[nodiscard]] static inline bool isMainThread() noexcept;
    private:
        ThreadPool(int threadCount);
        /* Operations */
        void workerMainLoop(int thread_id) noexcept;
        /* Static Members */
        [[nodiscard]] static ThreadInfo& getThreadInfo() noexcept;
        [[nodiscard]] static int getNumProcesser() noexcept;
        [[nodiscard]] static int makeNumThread() noexcept;
    };

    inline int ThreadPool::getThreadID() noexcept {
        return getThreadInfo().id;
    }

    inline bool ThreadPool::isMainThread() noexcept {
        return getThreadID() == MainThreadID;
    }
}
