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
#pragma once

#include <memory>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <condition_variable>
#include <future>
#include <sys/sysinfo.h>
#include <Physica/Utils/Container/Array/Array.h>
#include "PackagedTask.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Eigen; https://eigen.tuxfamily.org/
     */
    class PHYSICA_API ThreadPool final {
    public:
        struct ThreadData {
            std::unique_ptr<std::thread> thread;
            std::queue<std::unique_ptr<Task>> queue;
            std::mutex queueMutex;

            ThreadData() : thread(), queue(), queueMutex() {}
        };

        struct ThreadInfo {
            int id;
            uint64_t numScheduled;
            uint64_t randState;
        };

        static int numThreadRequired;
    private:
        constexpr static size_t MainThreadID = std::numeric_limits<decltype(numThreadRequired)>::max();
        thread_local static std::unique_ptr<ThreadInfo> info;

        Utils::Array<ThreadData> thread_data;
        std::mutex poolMutex;
        std::condition_variable cond;
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
        std::future<typename std::invoke_result<Function, Args&&...>::type> schedule(Function func, Args&&... args) noexcept;
        std::unique_ptr<Task> steal();
        void waitExit();
        void restart();
        /* Getters */
        [[nodiscard]] int getNumThreads() const noexcept { return thread_data.getLength(); }
        /* Setters */
        void shouldExit() noexcept;
        /* Static Members */
        [[nodiscard]] static ThreadInfo& getThreadInfo() noexcept;
        [[nodiscard]] static inline bool isMainThread() noexcept;
        [[nodiscard]] static ThreadPool& getInstance() noexcept;
    private:
        ThreadPool(int threadCount);
        /* Operations */
        void workerMainLoop(int thread_id) noexcept;
        /* Static Members */
        [[nodiscard]] static inline int getNumProcesser() noexcept { return get_nprocs(); }
        [[nodiscard]] static inline int makeNumThread() noexcept;
        [[nodiscard]] static inline unsigned int threadRand(uint64_t& state) noexcept;
    };

    template<class Function, class... Args>
    std::future<typename std::invoke_result<Function, Args&&...>::type> ThreadPool::schedule(Function func, Args&&... args) noexcept {
        using ResultType = typename std::invoke_result<Function, Args&&...>::type;
        int schedule_to;
        auto& info = getThreadInfo();
        if (isMainThread())
            schedule_to = info.numScheduled % getNumThreads();
        else
            schedule_to = info.id;
        info.numScheduled += 1;

        std::packaged_task<ResultType()> task(std::bind(func, std::forward<Args>(args)...));
        auto result = task.get_future();
        auto& data = thread_data[schedule_to];
        data.queueMutex.lock();
        data.queue.emplace(new PackagedTask(std::move(task)));
        data.queueMutex.unlock();
        cond.notify_one();
        return result;
    }

    inline bool ThreadPool::isMainThread() noexcept {
        return getThreadInfo().id == MainThreadID;
    }

    inline int ThreadPool::makeNumThread() noexcept {
        const auto numProcesser = getNumProcesser();
        if (numThreadRequired == 0 || numThreadRequired > numProcesser)
            return numProcesser * 3 / 4;
        else
            return numThreadRequired;
    }

    inline unsigned int ThreadPool::threadRand(uint64_t& state) noexcept {
        uint64_t current = state;
        state = current * 6364136223846793005ULL + 0xda3e39cb94b95bdbULL;
        // Generate the random output (using the PCG-XSH-RS scheme)
        return static_cast<unsigned int>((current ^ (current >> 22U)) >> (22U + (current >> 61U)));
    }
}
