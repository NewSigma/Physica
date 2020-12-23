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
#include "Physica/Logger/LoggerRuntime.h"
#include "Physica/Logger/Logger/StdLogger.h"

namespace Physica::Logger {
    LoggerRuntime::LoggerRuntime()
            : buffer(1U << 20U)
            , logThread(&LoggerRuntime::logThreadMain, this)
            , shouldExit(false) {
        registerLogger(std::unique_ptr<AbstractLogger>(new StdLogger(std::clog)));
        registerLogger(std::unique_ptr<AbstractLogger>(new StdLogger(std::cout)));
        registerLogger(std::unique_ptr<AbstractLogger>(new StdLogger(std::cerr)));
    }

    LoggerRuntime::~LoggerRuntime() {
        if(logThread.joinable())
            logThread.join();
        for (auto& logger : loggerList)
            delete logger;
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

    void LoggerRuntime::registerLogInfo(const LogInfo& info) {
        logInfos.push_back(info);
    }

    void LoggerRuntime::logThreadMain() {
        //Format [11:49:23] [Physica:12|Info]: This is a log.
        using namespace std::chrono_literals;
        while(!shouldExit || !buffer.isEmpty()) {
            while(buffer.isEmpty()) {
                if(shouldExit)
                    return;
                std::this_thread::sleep_for(1s);
            }
            size_t loggerID;
            buffer.read(&loggerID);
            loggerList[loggerID]->log();
        }
    }
}