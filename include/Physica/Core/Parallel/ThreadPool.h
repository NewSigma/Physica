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
#pragma once

#ifdef __linux__
    #include <sys/sysinfo.h>
#else
    #include <windows.h>
#endif
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
        class ThreadQueue {
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

        struct Awaiter : public suspend_always {
            static void await_suspend(Handle) noexcept;
        };
    public:
        static int numThreadRequired;
        // Larger than any thread ID, main thread is not maintained by thread pool
        constexpr static int MainThreadID = std::numeric_limits<decltype(numThreadRequired)>::max();
    private:
        Array<ThreadQueue> queues;
        std::mutex poolMutex;
        std::condition_variable cond;
        bool exit = false;
    public:
        ThreadPool(const This&) = delete;
        ThreadPool(This&&) noexcept = delete;
        ~ThreadPool();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] Awaiter operator co_await() noexcept;
        /* Operations */
        void notify_one();
        void notify_all();
        [[nodiscard]] Handle steal() noexcept;

        void shouldExit() noexcept;
        void waitExit();
        void restart();
        /* Getters */
        [[nodiscard]] int getNumThreads() const noexcept { return (int)queues.getLength(); }
        /* Static Members */
        [[nodiscard]] static This& getInstance() noexcept;
        [[nodiscard]] static int getThreadID() noexcept;
        [[nodiscard]] static bool isMainThread() noexcept { return getThreadID() == MainThreadID; }
        static void spin() noexcept;
    private:
        ThreadPool(int numThreads);
        /* Operations */
        void workerMainLoop(int thread_id) noexcept;
    };
}
