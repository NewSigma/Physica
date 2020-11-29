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
#include <Physica/Utils/Cycler.h>
#include <sstream>
#include <iostream>
#include <cstring>
#include <cassert>
#include "Physica/Logger/LoggerRuntime.h"

namespace Physica::Logger {
    Logger* LoggerRuntime::stdoutLogger = nullptr;

    LoggerRuntime::LoggerRuntime()
            : buffer(1U << 20U)
            , logThread(&LoggerRuntime::logThreadMain, this)
            , shouldExit(false) {}

    LoggerRuntime::~LoggerRuntime() {
        if(logThread.joinable())
            logThread.join();
    }

    void LoggerRuntime::registerLogger(const LogInfo& info) {
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
            std::stringstream logString{};
            size_t logID;
            buffer.read(&logID);
            assert(logID > 0);
            const LogInfo& info = logInfos[logID - 1];
            /* Handle time */ {
                size_t timeStart;
                buffer.read(&timeStart);
                time_t now = std::time(nullptr);
                size_t timeEnd = Utils::Cycler::now();
                now -= Utils::Cycler::toSeconds(timeEnd - timeStart);
                auto localTime = std::localtime(&now);
                logString << '['
                          << localTime->tm_hour
                          << ':'
                          << localTime->tm_min
                          << ':'
                          << localTime->tm_sec;
            }
            //Handle file, line and severity.
            logString << "] ["
                    << info.file;
            logString << ':'
                    << info.line
                    << '|'
                    << info.level
                    << "]: ";
            //Handle format
            const char* const format = info.format;
            size_t pos = 0;
            while(format[pos] != '\0') {
                if(format[pos] == '%') {
                    switch(format[pos + 1]) {
                        case '%':
                            logString << '%';
                        case 'c':
                            char temp;
                            buffer.read(&temp);
                            logString << static_cast<char>(temp);
                            break;
                        case 's':
                            //Bug: the string must be allocated statically.
                            char* str;
                            buffer.read(&str);
                            logString << str;
                            break;
                        case 'd':
                        case 'i':
                            int i;
                            buffer.read(&i);
                            logString << i;
                            break;
                        case 'o':
                        case 'x':
                        case 'X':
                        case 'u':
                        case 'f':
                        case 'F':
                        case 'a':
                        case 'A':
                        case 'g':
                        case 'G':
                            printf("[%s:%d|Fatal]: Logger not completely implemented.", __FILENAME__, __LINE__);
                            exit(EXIT_FAILURE);
                        case 'p':
                            void* p;
                            buffer.read(&p);
                            logString << p;
                            break;
                        default:
                            //The format should have been checked be format analyzer.
                            printf("[%s:%d|Fatal]: The code should not execute this sentence, please check your code."
                                   , __FILENAME__, __LINE__);
                            exit(EXIT_FAILURE);
                    }
                    ++pos;
                }
                else {
                    logString << format[pos];
                }
                ++pos;
            }
            //Output
            std::cout << logString.str() << std::endl;
        }
    }
}