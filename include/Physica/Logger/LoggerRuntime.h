/*
 * Copyright 2020 WeiBo He.
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
#ifndef PHYSICA_LOGGERRUNTIME_H
#define PHYSICA_LOGGERRUNTIME_H

#include <vector>
#include <thread>
#include "Physica/Utils/RingBuffer.h"
#include "LoggerType.h"
#include "Logger.h"

namespace Physica::Logger {
    class LoggerRuntime {
    public:
        constexpr static const char* __restrict levelString[4] = { "Fatal", "Warning", "Info", "Debug" };
    private:
        Utils::RingBuffer buffer;
        /*!
         * Store info of logged logs.
         */
        std::vector<LogInfo> logInfos;
        /*!
         * The thread used to output logs.
         */
        std::thread logThread;
        /*!
         * Does log thread should exit.
         */
        bool shouldExit;
    public:
        LoggerRuntime(const LoggerRuntime&) = delete;
        LoggerRuntime(LoggerRuntime&&) noexcept = delete;
        ~LoggerRuntime();
        /* Operators */
        LoggerRuntime& operator=(const LoggerRuntime&) = delete;
        LoggerRuntime& operator=(LoggerRuntime&&) noexcept = delete;
        /* Operations */
        void registerLogger(const LogInfo& info);
        void loggerShouldExit() { shouldExit = true; }
        /* Getters */
        [[nodiscard]] bool isIDRegistered(size_t id) const noexcept { return id != 0 && logInfos.size() >= id; }
        [[nodiscard]] size_t getNextLogID() const noexcept { return logInfos.size(); }
        [[nodiscard]] Utils::RingBuffer& getBuffer() noexcept { return buffer; }
        /* Static Members */
        static inline LoggerRuntime& getInstance();
    private:
        LoggerRuntime();

        void logThreadMain();
    };

    inline LoggerRuntime& LoggerRuntime::getInstance() {
        static LoggerRuntime runtime{};
        return runtime;
    }

    template<typename T1, typename... Ts>
    inline void writeArgs(T1 head, Ts... args);

    inline void writeArgs();

    template<typename... Ts>
    void log(size_t logID, Ts... args);
}

#include "LoggerRuntimeImpl.h"

#endif
