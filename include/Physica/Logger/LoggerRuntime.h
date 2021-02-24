/*
 * Copyright 2020-2021 WeiBo He.
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

#include <cassert>
#include <vector>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <memory>
#include <cstring>
#include "LogBuffer.h"
#include "LoggerType.h"
#include "Physica/Logger/Logger/AbstractLogger.h"
#include "LoggerTimer.h"

namespace Physica::Logger {
    /**
     * Three loggers will be created after initialized: std::clog, std::cout, std::cerr,
     * whose ids are 0, 1 and 2.
     *
     * The ids are arranged to equal to the file descriptors of stdio.
     */
    class LoggerRuntime {
    public:
        constexpr static const char* __restrict levelString[4] = { "Fatal", "Warning", "Info", "Debug" };
        constexpr static size_t unassignedLogID = 0;
    private:
        thread_local static LogBuffer* threadLogBuffer;

        LoggerTimer timer;
        std::vector<LogBuffer*> bufferList;
        std::mutex bufferListMutex;
        /**
         * ID of buffer being logged.
         */
        int processingBufferID;
        /**
         * Store all registered loggers.
         */
        std::vector<AbstractLogger*> loggerList;
        /**
         * Store info of logged logs.
         */
        std::vector<LogInfo> logInfos;
        /**
         * The thread used to output logs.
         */
        std::thread logThread;
        /**
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
        size_t registerLogger(std::unique_ptr<AbstractLogger>&& logger);
        void registerLogInfo(const LogInfo& info) { logInfos.push_back(info); }
        void loggerShouldExit() { shouldExit = true; }
        /* Getters */
        [[nodiscard]] const LoggerTimer& getTimer() const noexcept { return timer; }
        [[nodiscard]] AbstractLogger& getLogger(size_t index) const { assert(index < loggerList.size()); return *loggerList[index]; }
        [[nodiscard]] size_t getNextLogID() const noexcept { return logInfos.size(); }
        [[nodiscard]] LogInfo getLogInfo(size_t index) const { return logInfos[index]; }
        [[nodiscard]] Utils::RingBuffer& getBuffer();
        /* Setters */
        void releaseBuffer() noexcept { threadLogBuffer->schedualDelete(); threadLogBuffer = nullptr; }
        /* Static Members */
        static inline LoggerRuntime& getInstance();
    private:
        LoggerRuntime();

        void logThreadMain();
        int getNextBufferToLog() noexcept;
    };

    inline AbstractLogger& getLogger(size_t index) {
        return LoggerRuntime::getInstance().getLogger(index);
    }

    inline LoggerRuntime& LoggerRuntime::getInstance() {
        static LoggerRuntime runtime{};
        return runtime;
    }

    template<typename... Ts>
    void log(const ArgType* p_args, size_t logID, Ts... args);
}

#include "LoggerRuntimeImpl.h"

#define Log(loggerID, severity, format, ...)                                                        \
    do {                                                                                            \
        using namespace Physica::Logger;                                                            \
        AbstractLogger& logger = LoggerRuntime::getInstance().getLogger(loggerID);                  \
        LogLevel level = logger.getCurrentLevel();                                                  \
        if(level >= LogLevel::severity) {                                                           \
            if(false)                                                                               \
                Physica::Logger::checkFormat(format, ##__VA_ARGS__);                                \
                                                                                                    \
            constexpr size_t argCount = getArgCount(format);                                        \
            static constexpr std::array<ArgType, argCount> argArray                                 \
                                                        = analyzeFormatString<argCount>(format);    \
            static size_t logID = LoggerRuntime::unassignedLogID;                                   \
                                                                                                    \
            if(logID == LoggerRuntime::unassignedLogID) {                                           \
                constexpr LogInfo info{                                                             \
                        LoggerRuntime::levelString[static_cast<int>(LogLevel::severity)],           \
                        format,                                                                     \
                        __FILENAME__,                                                               \
                        __LINE__,                                                                   \
                        argArray.data(),                                                            \
                        argCount};                                                                  \
                LoggerRuntime::getInstance().registerLogInfo(info);                                 \
                logID = LoggerRuntime::getInstance().getNextLogID();                                \
            }                                                                                       \
            writeArgs(argArray.begin(), static_cast<size_t>(loggerID));                             \
            Physica::Logger::log(argArray.begin(), logID, ##__VA_ARGS__);                           \
        }                                                                                           \
    } while(false)

#ifndef NDEBUG
    #define Debug(loggerID, format, ...) Log(loggerID, Debug, format, ##__VA_ARGS__)
#else
    #define Debug(loggerID, format, ...) do {} while(false)
#endif

#define Info(loggerID, format, ...) Log(loggerID, Info, format, ##__VA_ARGS__)

#define Warning(loggerID, format, ...) Log(loggerID, Warning, format, ##__VA_ARGS__)

#define Fatal(loggerID, format, ...) do { Log(loggerID, Fatal, format, ##__VA_ARGS__); exit(EXIT_FAILURE); } while(false)
/**
 * Use Error instead of Fatal when a system call fails.
 */
#define Error(loggerID, message) Fatal(loggerID, "%s: %s", message, strerror(errno))
