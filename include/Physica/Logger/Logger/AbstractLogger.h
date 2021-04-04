/*
 * Copyright 2020-2021 WeiBo He.
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
#pragma once

#include <cstdarg>
#include <iostream>
#include <array>
#include "Physica/Config.h"
#include "Physica/Logger/LogBuffer.h"
#include "Physica/Logger/LoggerType.h"
#include "Physica/Logger/FormatAnalyzer.h"

namespace Physica::Logger {
    /*!
     * Abstract father class for all loggers.
     */
    class AbstractLogger {
        static LogLevel globalLevel;
    public:
        LogLevel localLevel;
    public:
        explicit AbstractLogger(LogLevel level = LogLevel::Global) : localLevel(level) {}
        AbstractLogger(const AbstractLogger& logger) = delete;
        AbstractLogger(AbstractLogger&& logger) noexcept = delete;
        virtual ~AbstractLogger() = default;
        /* Operators */
        AbstractLogger& operator=(const AbstractLogger&) = delete;
        AbstractLogger& operator=(AbstractLogger&&) noexcept = delete;
        /* Operations */
        virtual void log(LogBuffer& buffer) = 0;
        /* Getters */
        [[nodiscard]] inline LogLevel getCurrentLevel() const;
        /* Static members */
        [[nodiscard]] static LogLevel getGlobalLevel() noexcept { return globalLevel; }
        static inline void setGlobalLevel(LogLevel level) noexcept;
    };

    inline LogLevel AbstractLogger::getCurrentLevel() const {
        return localLevel == LogLevel::Global ? globalLevel : localLevel;
    }

    inline void AbstractLogger::setGlobalLevel(LogLevel level) noexcept {
        globalLevel = level == LogLevel::Global ? LogLevel::Off : level;
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
