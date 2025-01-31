/*
 * Copyright 2020-2025 Weibo He.
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
#include "Physica/Logger/LoggerRuntime.h"

namespace Physica::Logger {
    /*!
     * Abstract father class for all loggers.
     */
    class AbstractLogger {
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
        [[nodiscard]] inline LogLevel getCurrentLevel() const noexcept;
    };

    inline LogLevel AbstractLogger::getCurrentLevel() const noexcept {
        return localLevel == LogLevel::Global ? LoggerRuntime::globalLevel : localLevel;
    }
}
