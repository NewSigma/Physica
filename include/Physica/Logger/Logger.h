/*
 * Copyright 2026 Weibo He.
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
#pragma once // FIXME: Make it a header unit once cmake support it. Ref: https://gitlab.kitware.com/cmake/cmake/-/issues/25293

#include <array>

import Physica.Logger;

#define Log(loggerID, severity, format, ...)                                                        \
    do {                                                                                            \
        using namespace Physica;                                                                    \
        auto& logger = LoggerRuntime::getInstance().getLogger(loggerID);                            \
        LogLevel level = logger.getCurrentLevel();                                                  \
        if(level >= LogLevel::severity) {                                                           \
            constexpr size_t argCount = [](auto... args) consteval {                                \
                return sizeof...(args);                                                             \
            }(__VA_ARGS__);                                                                       \
            constexpr static auto args = analyzeFormatString<argCount>(format);                     \
            static size_t logID = LoggerRuntime::unassignedLogID;                                   \
                                                                                                    \
            if(logID == LoggerRuntime::unassignedLogID) {                                           \
                constexpr LogInfo info{                                                             \
                        LoggerRuntime::levelString[static_cast<int>(LogLevel::severity) - 1],       \
                        format,                                                                     \
                        getFileName(__FILE__),                                                      \
                        __LINE__,                                                                   \
                        args.data(),                                                                \
                        args.size()};                                                               \
                LoggerRuntime::getInstance().registerLogInfo(info);                                 \
                logID = LoggerRuntime::getInstance().getNextLogID();                                \
            }                                                                                       \
            writeArgs(args.begin(), static_cast<size_t>(loggerID));                                 \
            Physica::log(args.begin(), logID, ##__VA_ARGS__);                                       \
        }                                                                                           \
    } while(false)

#ifndef NDEBUG
    #define Debug(loggerID, format, ...) Log(loggerID, Debug, format, ##__VA_ARGS__)
#else
    #define Debug(loggerID, format, ...) do {} while(false)
#endif

#define Info(loggerID, format, ...) Log(loggerID, Info, format, ##__VA_ARGS__)

#define Warn(loggerID, format, ...) Log(loggerID, Warn, format, ##__VA_ARGS__)

#define Fatal(loggerID, format, ...)                    \
    do {                                                \
        using namespace Physica;                        \
        Log(loggerID, Fatal, format, ##__VA_ARGS__);    \
        std::abort();                                   \
    } while(false)
/**
 * Use Error instead of Fatal when a system call fails.
 */
#define Error(loggerID, message) Fatal(loggerID, "%s: %s", message, strerror(errno))

// TODO: A better logger for both C++ and CUDA is required.
#define cuDebug(x) do { printf("[] [Debug] [%s:%d]: %s\n", __FILENAME__, __LINE__, x); } while(false)
#define cuWarning(x) do { printf("[] [Warning] [%s:%d]: %s\n", __FILENAME__, __LINE__, x); } while(false)
#define cuCritical(x) do { printf("[] [Critical] [%s:%d]: %s\n", __FILENAME__, __LINE__, x); } while(false)
#define cuFatal(x)                                                    \
    do {                                                              \
        printf("[] [Fatal] [%s:%d]: %s\n", __FILENAME__, __LINE__, x);\
        cudaDeviceReset();                                            \
        exit(EXIT_FAILURE);                                           \
    } while(false)
#define cuInfo(x) do { printf("[] [Info] [%s:%d]: %s\n", __FILENAME__, __LINE__, x); } while(false)
