/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#ifdef __linux__
    #include <sys/sysinfo.h>
#else
    #include <windows.h>
#endif
#include <mutex>
#include <condition_variable>
#include <limits>
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    /**
     * A task-opaque thread pool that supports stealing-based load balancing.
     *
     * TODO: Topology sensing
     *
     * Reference:
     * [1] Eigen; https://eigen.tuxfamily.org
     */
    class PHYSICA_API ThreadPool final {
        using This = ThreadPool;
        using Handle = std::coroutine_handle<>;

        class ThreadQueue;

        // TODO: Align with the split of std::execution once we bump to CXX26
        struct ImplicitScheduler : public suspend_always {
            static void await_suspend(Handle) noexcept;
            static void on_wait(Handle) noexcept;
        };

        struct Scheduler : public ImplicitScheduler {
            using Base = ImplicitScheduler;

            static Base implicit() noexcept;
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
        [[nodiscard]] Scheduler operator co_await() noexcept;
        /* Operations */
        void notify_one() noexcept;
        void notify_all() noexcept;
        [[nodiscard]] Handle steal() noexcept;

        void resize(int numThreads);
        void waitExit();
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
        void spawnWorkers();
        void workerMainLoop(int thread_id) noexcept;
    };
}
