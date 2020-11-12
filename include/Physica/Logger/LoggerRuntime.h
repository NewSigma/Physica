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
#ifndef PHYSICA_LOGGERRUNTIME_H
#define PHYSICA_LOGGERRUNTIME_H

#include "LoggerType.h"

namespace Physica::Logger {
    class LoggerRuntime {
    public:
        constexpr static const char* __restrict levelString[4] = { "Fatal", "Warning", "Info", "Debug" };
        static LogLevel globalLevel;
    private:
        /*!
         * The buffer between the producer and the consumer.
         */
        char buffer[1U << 20U];
    public:
        /* Static Members */
        static inline LoggerRuntime& getInstance();
        [[nodiscard]] static LogLevel getGlobalLevel() { return globalLevel; }
        static inline void setGlobalLevel(LogLevel level);
    };

    inline LoggerRuntime& LoggerRuntime::getInstance() {
        static LoggerRuntime runtime{};
        return runtime;
    }

    inline void LoggerRuntime::setGlobalLevel(LogLevel level) {
        globalLevel = level == LogLevel::Global ? LogLevel::Off : level;
    }
}

#endif
