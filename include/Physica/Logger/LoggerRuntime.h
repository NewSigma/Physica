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
#include "AbstractLogger.h"

namespace Physica::Logger {
    class LoggerRuntime {
    public:
        constexpr static const char* __restrict levelString[4] = { "Fatal", "Warning", "Info", "Debug" };
        constexpr static const size_t unassignedLogID = 0;
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
        [[nodiscard]] size_t getNextLogID() const noexcept { return logInfos.size(); }
        [[nodiscard]] Utils::RingBuffer& getBuffer() noexcept { return buffer; }
        /* Static Members */
        static inline AbstractLogger& getStdoutLogger();
        static inline LoggerRuntime& getInstance();
    private:
        LoggerRuntime();

        void logThreadMain();
    };

    inline AbstractLogger& LoggerRuntime::getStdoutLogger() {
        static AbstractLogger stdoutLogger{};
        return stdoutLogger;
    }

    inline LoggerRuntime& LoggerRuntime::getInstance() {
        static LoggerRuntime runtime{};
        return runtime;
    }
    /*!
     * It is same to call LoggerRuntime::getStdoutLogger(), but this function is more convenient.
     */
    inline AbstractLogger& getStdoutLogger() { return LoggerRuntime::getStdoutLogger(); }

    template<typename... Ts>
    void log(const ArgType* p_args, size_t logID, Ts... args);
}

#define Log(logger, severity, format, ...)                                                          \
    do {                                                                                            \
        using namespace Physica::Logger;                                                            \
        LogLevel level = logger.getCurrentLevel();                                                  \
        if(level >= LogLevel::severity) {                                                           \
            if(false)                                                                               \
                Physica::Logger::checkFormat(format, ##__VA_ARGS__);                                \
            static size_t logID = LoggerRuntime::unassignedLogID;                                   \
                                                                                                    \
            constexpr size_t argCount = getArgCount(format);                                        \
            static constexpr std::array<ArgType, argCount> argArray                                 \
                                                        = analyzeFormatString<argCount>(format);    \
            if(logID == LoggerRuntime::unassignedLogID) {                                           \
                constexpr LogInfo info{                                                             \
                        LoggerRuntime::levelString[static_cast<int>(LogLevel::severity)],           \
                        format,                                                                     \
                        __FILENAME__,                                                               \
                        __LINE__,                                                                   \
                        argArray.data(),                                                            \
                        argCount};                                                                  \
                LoggerRuntime::getInstance().registerLogger(info);                                  \
                logID = LoggerRuntime::getInstance().getNextLogID();                                \
            }                                                                                       \
            Physica::Logger::log(argArray.begin(), logID, ##__VA_ARGS__);                           \
        }                                                                                           \
    } while(false)

#ifndef NDEBUG
    #define Debug(logger, format, ...) Log(logger, Debug, format, ##__VA_ARGS__)
#else
    #define Debug(logger, format, ...) do {} while(false)
#endif

#define Info(logger, format, ...) Log(logger, Info, format, ##__VA_ARGS__)

#define Warning(logger, format, ...) Log(logger, Warning, format, ##__VA_ARGS__)

#define Fatal(logger, format, ...) do { Log(logger, Fatal, format, ##__VA_ARGS__); exit(EXIT_FAILURE); } while(false)

#include "LoggerRuntimeImpl.h"

#endif
