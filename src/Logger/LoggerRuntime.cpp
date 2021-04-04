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
#include <iostream>
#include <cstring>
#include <unistd.h>
#include "Physica/Logger/LoggerRuntime.h"
#include "Physica/Logger/Logger/StdLogger.h"

namespace Physica::Logger {
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
        getBuffer();

        logThread = std::thread(&LoggerRuntime::logThreadMain, this);
    }

    LoggerRuntime::~LoggerRuntime() {
        if (logThread.joinable())
            logThread.join();
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

    Utils::RingBuffer& LoggerRuntime::getBuffer() {
        if (shouldExit)
            Fatal(STDERR_FILENO, "Try to append log to closed LoggerRuntime.");
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
                getNextBufferToLog();
            }

            LogBuffer& buffer = *bufferList[processingBufferID];
            while (!buffer.isEmpty()) {
                size_t loggerID;
                buffer.read(&loggerID);
                loggerList[loggerID]->log(buffer);
            }
            getNextBufferToLog();
        }
    }
    /**
     * Return true if all buffers are empty.
     */
    int LoggerRuntime::getNextBufferToLog() noexcept {
        if (processingBufferID == -1)
            return 0;
        size_t size = bufferList.size();
        int i = static_cast<int>((processingBufferID + 1) % size);
        for (; i != processingBufferID; i = (processingBufferID + 1) % size) {
            LogBuffer* buffer = bufferList[i];
            if (buffer->isEmpty()) {
                if (buffer->getShouldDelete()) {
                    std::unique_lock<std::mutex> lock(bufferListMutex);
                    delete buffer;
                    bufferList.erase(bufferList.begin() + i);
                    --size;
                }
            }
            else {
                processingBufferID = i;
                return processingBufferID;
            }
        }
        processingBufferID = -1;
        return processingBufferID;
    }
}