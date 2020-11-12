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
#ifndef PHYSICA_LOGGER_H
#define PHYSICA_LOGGER_H

#include <cstdarg>
#include <iostream>
#include "Physica/Config.h"
#include "LoggerType.h"
#include "LoggerRuntime.h"
#include "FormatAnalyzer.h"

namespace Physica::Logger {
    /*!
     * Abstract father class for all loggers.
     */
    class Logger {
    public:
        LogLevel localLevel;
    public:
        explicit Logger(LogLevel level = LogLevel::Global) : localLevel(level) {}
        Logger(const Logger& logger) = default;
        Logger(Logger&& logger) noexcept = default;
        ~Logger() = default;
        /* Operators */
        Logger& operator=(const Logger&) = default;
        Logger& operator=(Logger&&) noexcept = default;
        /* Operations */
        static void log(const LogInfo* info, ...);
        /* Getters */
        [[nodiscard]] inline LogLevel getCurrentLevel() const;
    };

    inline LogLevel Logger::getCurrentLevel() const {
        return localLevel == LogLevel::Global ? LoggerRuntime::getGlobalLevel() : localLevel;
    }
    /*!
     * No-Op function that triggers the GNU preprocessor's format checker for
     * printf format strings and argument parameters.
     *
     * \param format
     *      printf format string
     * \param ...
     *      format parameters
     */
    inline void __attribute__ ((format (printf, 1, 2))) checkFormat(const char*, ...) {}
}

#define Log(severity, format, ...)                                                      \
    do {                                                                                \
        using namespace Physica::Logger;                                                \
        LogLevel level = logger.getCurrentLevel();                                      \
        if(level >= LogLevel::severity) {                                               \
            if(false)                                                                   \
                Physica::Logger::checkFormat(format, ##__VA_ARGS__);                    \
            constexpr LogInfo info{                                                     \
                    LoggerRuntime::levelString[static_cast<int>(LogLevel::severity)],   \
                    format,                                                             \
                    __FILENAME__,                                                       \
                    __LINE__};                                                          \
            logger.log(&info, ##__VA_ARGS__);                                           \
        }                                                                               \
    } while(false)

#define Debug(logger, format, ...) Log(Debug, format, ##__VA_ARGS__)

#define Info(logger, format, ...) Log(Info, format, ##__VA_ARGS__)

#define Warning(logger, format, ...) Log(Warning, format, ##__VA_ARGS__)

#define Fatal(logger, format, ...) Log(Fatal, format, ##__VA_ARGS__)

#endif
