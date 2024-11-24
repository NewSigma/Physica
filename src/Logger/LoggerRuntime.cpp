/*
 * Copyright 2020-2024 Weibo He.
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
#include <iostream>
#include <cstring>
#include <csignal>
#include <unistd.h>
#include "Physica/Logger/LoggerRuntime.h"
#include "Physica/Logger/Logger/StdLogger.h"

namespace Physica::Logger {
    LogLevel LoggerRuntime::globalLevel = LogLevel::Info;
    thread_local LogBuffer* LoggerRuntime::threadLogBuffer = nullptr;

    LoggerRuntime::LoggerRuntime()
            : timer()
            , bufferList()
            , bufferListMutex()
            , processingBufferID(0)
            , shouldExit(false) {
        registerLogger(std::unique_ptr<AbstractLogger>(new StdLogger(std::clog)));
        registerLogger(std::unique_ptr<AbstractLogger>(new StdLogger(std::cout)));
        registerLogger(std::unique_ptr<AbstractLogger>(new StdLogger(std::cerr)));
        //Init buffer for current thread or logThread will try to access a empty bufferList
        std::ignore = getBuffer();

        logThread = std::thread(&LoggerRuntime::logThreadMain, this);

        if (std::signal(SIGABRT, abort_handler) == SIG_ERR) {
            Debug(STDERR_FILENO, "Implementation forbid us from handling SIGABRT");
        }
    }

    LoggerRuntime::~LoggerRuntime() {
        waitExit();
        for (auto& logger : loggerList)
            delete logger;
        for (LogBuffer* buffer : bufferList)
            delete buffer;
    }
    /**
     * \param logger
     * Pointer to a logger thar is allocated on heap.
     *
     * \return
     * The id of the registered logger.
     */
    size_t LoggerRuntime::registerLogger(std::unique_ptr<AbstractLogger>&& logger) {
        auto nextID = loggerList.size();
        loggerList.push_back(logger.release());
        return nextID;
    }

    void LoggerRuntime::waitExit() noexcept {
        loggerShouldExit();
        if (logThread.joinable())
            logThread.join();
    }

    Core::RingBuffer& LoggerRuntime::getBuffer() {
        if (shouldExit) [[unlikely]]
            throw std::runtime_error("[Error]: Try to append log to closed LoggerRuntime");
        if (threadLogBuffer == nullptr) {
            threadLogBuffer = new LogBuffer(DefaultBufferSize);
            std::unique_lock<std::mutex> lock(bufferListMutex);
            bufferList.push_back(threadLogBuffer);
        }
        return *threadLogBuffer;
    }

    void LoggerRuntime::logThreadMain() {
        //Format [11:49:23] [Physica:12|Info]: This is a log.
        using namespace std::chrono_literals;
        
        while (!shouldExit || (processingBufferID >= 0)) {
            while (processingBufferID < 0) {
                if(shouldExit)
                    return;
                std::this_thread::sleep_for(1s);
                findNextBufferToLog();
            }

            bufferListMutex.lock();
            LogBuffer& buffer = *bufferList[processingBufferID];
            bufferListMutex.unlock();

            while (!buffer.isEmpty()) {
                size_t loggerID;
                buffer.read(&loggerID);
                loggerList[loggerID]->log(buffer);
            }
            findNextBufferToLog();
        }
    }

    void LoggerRuntime::findNextBufferToLog() noexcept {
        std::unique_lock<std::mutex> lock(bufferListMutex);
        std::vector<LogBuffer*> newBufferList{};
        newBufferList.reserve(bufferList.size());
        for (auto ite = bufferList.begin(); ite != bufferList.end(); ++ite) {
            auto buffer = *ite;
            if (!(buffer->isEmpty() && buffer->getShouldDelete()))
                newBufferList.push_back(buffer);
        }
        bufferList.swap(newBufferList);

        const size_t size = bufferList.size();
        for (size_t i = 0; i < size; ++i) {
            processingBufferID = static_cast<int>((processingBufferID + 1) % size);
            auto* buffer = bufferList[processingBufferID];
            if (!buffer->isEmpty())
                return;
        }
        processingBufferID = -1;
    }

    void LoggerRuntime::abort_handler(int) noexcept {
        getInstance().~LoggerRuntime();
        std::_Exit(EXIT_FAILURE);
    }
}
