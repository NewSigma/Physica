/*
 * Copyright 2021 WeiBo He.
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
#include <sstream>
#include <iostream>
#include <cassert>
#include <ctime>
#include "Physica/Logger/LogBuffer.h"
#include "Physica/Logger/LoggerRuntime.h"
#include "Physica/Utils/Cycler.h"

namespace Physica::Logger {
    std::string LogBuffer::makeMsgString() {
        std::stringstream logString{};
        size_t logID;
        LoggerRuntime& runtime = LoggerRuntime::getInstance();
        read(&logID);
        assert(logID > 0);
        const LogInfo& info = runtime.getLogInfo(logID - 1);
        /* Handle time */ {
            const LoggerTimer& timer = runtime.getTimer();
            uint64_t cycle;
            read(&cycle);
            time_t now = timer.toTime(cycle);
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
        return logString.str() + formatToString(info.format);
    }

    std::string LogBuffer::formatToString(const char* __restrict format) {
        std::stringstream logString{};
        size_t pos = 0;
        while(format[pos] != '\0') {
            if(format[pos] == '%') {
                switch(format[pos + 1]) {
                    case '%':
                        logString << '%';
                        break;
                    case 'c':
                        char temp;
                        read(&temp);
                        logString << static_cast<char>(temp);
                        break;
                    case 's': {
                        size_t strLength;
                        read(&strLength);
                        char* str = new char[strLength + 1];
                        char* p = str;
                        for(size_t i = 0; i < strLength; ++i)
                            read(p++);
                        str[strLength] = '\0';
                        logString << str;
                        break;
                    }
                    case 'd':
                    case 'i':
                        int i;
                        read(&i);
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
                        read(&p);
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
        return logString.str();
    }

    size_t LogBuffer::getMsgSize(const char* __restrict format) const {
        size_t pos = 0;
        size_t result = 0;
        while(format[pos] != '\0') {
            if(format[pos] == '%') {
                switch(format[pos + 1]) {
                    case '%':
                        break;
                    case 'c':
                        result += sizeof(char);
                        break;
                    case 's': {
                        size_t strLength;
                        cread(&strLength, result);
                        result += sizeof(size_t) + strLength;
                        break;
                    }
                    case 'd':
                    case 'i':
                        result += sizeof(int);
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
                        result += sizeof(void*);
                        break;
                    default:
                        //The format should have been checked be format analyzer.
                        printf("[%s:%d|Fatal]: The code should not execute this sentence, please check your code."
                                , __FILENAME__, __LINE__);
                        exit(EXIT_FAILURE);
                }
                ++pos;
            }
            ++pos;
        }
        return result;
    }
}
