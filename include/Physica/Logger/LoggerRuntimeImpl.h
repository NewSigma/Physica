/*
 * Copyright 2020-2024 Weibo He.
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

#include <cstring>
#include "Physica/Core/Utils/Cycler.h"

namespace Physica {
    template<typename T1, typename... Ts>
    inline
    std::enable_if<!std::is_same<T1, char*>::value
                            && !std::is_same<T1, const char*>::value>::type
    writeArgs(const ArgType* p_args, T1 head, Ts... args) {
        LoggerRuntime::getInstance().getBuffer().write(head);
        writeArgs(p_args + 1, args...);
    }

    template<typename T1, typename... Ts>
    inline
    std::enable_if<std::is_same<T1, char*>::value
                            || std::is_same<T1, const char*>::value>::type
    writeArgs(const ArgType* p_args, T1 head, Ts... args) {
        RingBuffer& buffer = LoggerRuntime::getInstance().getBuffer();
        if(*p_args == ArgType::s) {
            size_t strLength = std::strlen(head);
            buffer.write(strLength);
            buffer.writeBytes(head, strLength);
        }
        else {
            buffer.write(head);
        }
        writeArgs(p_args + 1, args...);
    }

    inline void writeArgs(const ArgType* p_args) { (void)p_args; } //Do nothing
    /**
     * \param p_args
     * Pointer to the first element of the ArgType array.
     *
     * \param logID
     * The id of the log to be logged.
     *
     * \param args
     * Arg pack of this log.
     */
    template<typename... Ts>
    void log(const ArgType* p_args, size_t logID, Ts... args) {
        assert(logID <= LoggerRuntime::getInstance().getNextLogID() && "[Error]: The log ID is not registered");
        size_t time = Cycler::now();
        writeArgs(p_args, logID);
        writeArgs(p_args, time);
        writeArgs(p_args, args...);
    }
}
