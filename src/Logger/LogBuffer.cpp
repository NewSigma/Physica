/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Logger/LogBuffer.h"
#include <cassert>
#include <ctime>
#include <iostream>
#include <sstream>
#include "Physica/Core/Utils/Cycler.h"
#include "Physica/Logger/LoggerRuntime.h"

using namespace Physica;

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
        const timeval t = timer.toTime(cycle);
        auto localTime = std::localtime(&t.tv_sec);
        logString << '['
                  << localTime->tm_hour
                  << ':'
                  << localTime->tm_min
                  << ':'
                  << localTime->tm_sec
                  << ':'
                  << (t.tv_usec / 1000);
    }
    // Handle file, line and severity.
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
    while (format[pos] != '\0') {
        if (format[pos] == '%') {
            switch (format[pos + 1]) {
            case '%':
                logString << '%';
                break;
            case 'c': {
                char temp;
                read(&temp);
                logString << static_cast<char>(temp);
                break;
            }
            case 's': {
                size_t strLength;
                read(&strLength);
                char* str = new char[strLength + 1];
                readBytes(str, strLength);
                str[strLength] = '\0';
                logString << str;
                delete[] str;
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
                printf("[%s:%d|Fatal]: Logger not completely implemented.", getFileName(__FILE__), __LINE__);
                exit(EXIT_FAILURE);
            case 'u': {
                unsigned int temp;
                read(&temp);
                logString << temp;
                break;
            }
            case 'U': {
                unsigned long temp;
                read(&temp);
                logString << temp;
                break;
            }
            case 'f':
                double f;
                read(&f);
                logString << f;
                break;
            case 'F':
            case 'a':
            case 'A':
            case 'g':
            case 'G':
                printf("[%s:%d|Fatal]: Logger not completely implemented.", getFileName(__FILE__), __LINE__);
                exit(EXIT_FAILURE);
            case 'p':
                void* p;
                read(&p);
                logString << p;
                break;
            default:
                // The format should have been checked by the format analyzer.
                printf("[%s:%d|Fatal]: The code should not execute this sentence, please check your code.", getFileName(__FILE__), __LINE__);
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
    while (format[pos] != '\0') {
        if (format[pos] == '%') {
            switch (format[pos + 1]) {
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
                printf("[%s:%d|Fatal]: Logger not completely implemented.", getFileName(__FILE__), __LINE__);
                exit(EXIT_FAILURE);
            case 'p':
                result += sizeof(void*);
                break;
            default:
                // The format should have been checked be format analyzer.
                printf("[%s:%d|Fatal]: The code should not execute this sentence, please check your code.", getFileName(__FILE__), __LINE__);
                exit(EXIT_FAILURE);
            }
            ++pos;
        }
        ++pos;
    }
    return result;
}
