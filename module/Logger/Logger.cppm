/*
 * Copyright 2020-2026 Weibo He.
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
module;

#include <vector>
#include <unordered_map>
#include <thread>
#include "Physica/Core/Utils/Container/RingBuffer.h"
#include "Physica/Core/Utils/Cycler.h"

export module Physica.Logger;

export import Physica.Logger.FormatAnalyzer;
export import Physica.Logger.LoggerType;
export import Physica.Logger.AbstractLogger;
import Physica.Logger.LogBuffer;
import Physica.Logger.LoggerTimer;

export namespace Physica {
    /**
     * Three loggers will be created after initialized: std::clog, std::cout, std::cerr,
     * whose ids are 0, 1 and 2.
     *
     * The ids are arranged to equal to the file descriptors of stdio.
     * 
     * The implementation is based on NanaLog[1]
     * 
     * References:
     * [1] Yang. Stephen,Park. Seo Jin,Ousterhout. John. NanoLog: A nanosecond scale logging system[J]. Proceedings of the 2018 USENIX Annual Technical Conference, USENIX ATC 2018:335-349
     */
    class PHYSICA_API LoggerRuntime {
    public:
        constexpr static const char* __restrict levelString[4] = { "Fatal", "Warn", "Info", "Debug" };
        constexpr static size_t unassignedLogID = 0;
        static LogLevel globalLevel;
    private:
        constexpr static size_t DefaultBufferSize = 1U << 20U;
        thread_local static LogBuffer* threadLogBuffer;

        LoggerTimer timer;
        std::vector<LogBuffer*> bufferList;
        std::mutex bufferListMutex;
        /**
         * ID of buffer being logged.
         */
        int processingBufferID;
        std::vector<AbstractLogger*> loggerList;
        std::vector<LogInfo> logInfos;
        /**
         * The thread used to output logs.
         */
        std::thread logThread;
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
        void loggerShouldExit() noexcept { shouldExit = true; }
        void waitExit() noexcept;
        /* Getters */
        [[nodiscard]] const LoggerTimer& getTimer() const noexcept { return timer; }
        [[nodiscard]] AbstractLogger& getLogger(size_t index) const { assert(index < loggerList.size()); return *loggerList[index]; }
        [[nodiscard]] size_t getNextLogID() const noexcept { return logInfos.size(); }
        [[nodiscard]] LogInfo getLogInfo(size_t index) const { return logInfos[index]; }
        [[nodiscard]] RingBuffer& getBuffer();
        /* Setters */
        void releaseBuffer() noexcept { threadLogBuffer->schedualDelete(); threadLogBuffer = nullptr; }
        /* Static Members */
        [[nodiscard]] static LoggerRuntime& getInstance() noexcept;
    private:
        LoggerRuntime();
        /* Operations */
        void logThreadMain();
        void findNextBufferToLog() noexcept;
        /* Static members */
        [[noreturn]] static void abort_handler(int sig) noexcept;
    };

    inline AbstractLogger& getLogger(size_t index) {
        return LoggerRuntime::getInstance().getLogger(index);
    }

    template<typename T1, typename... Ts>
    void writeArgs(const ArgType* p_args, T1 head, Ts... args);
    /**
     * \param p_args
     * Pointer to the first element of the ArgType array.
     */
    void log(const ArgType* p_args, size_t logID, auto... args);
}

namespace Physica {
    inline void writeArgs(const ArgType*) noexcept {}

    template<typename T1, typename... Ts>
    void writeArgs(const ArgType* p_args, T1 head, Ts... args) {
        constexpr bool isString = std::same_as<T1, char*> || std::same_as<T1, const char*>;
        RingBuffer& buffer = LoggerRuntime::getInstance().getBuffer();
        if constexpr (isString) {
            if (*p_args == ArgType::s) {
                size_t strLength = std::strlen(head);
                buffer.write(strLength);
                buffer.writeBytes(head, strLength);
            }
            else
                buffer.write(head);
        }
        else
            buffer.write(head);
        writeArgs(p_args + 1, args...);
    }
    /**
     * \param p_args
     * Pointer to the first element of the ArgType array.
     */
    void log(const ArgType* p_args, size_t logID, auto... args) {
        assert(logID <= LoggerRuntime::getInstance().getNextLogID() && "[Error]: The log ID is not registered");
        size_t time = Cycler::now();
        writeArgs(p_args, logID);
        writeArgs(p_args, time);
        writeArgs(p_args, args...);
    }
}
