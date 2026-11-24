/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "LoggerRuntime.h"

namespace Physica {
    inline void writeArgs(const ArgType*) noexcept {}

    template<typename T1, typename... Ts>
    void writeArgs(const ArgType* p_args, T1 head, Ts... args) {
        constexpr bool isString = std::same_as<T1, char*> || std::same_as<T1, const char*>;
        RingBuffer& buffer = LoggerRuntime::getInstance().getBuffer();
        if constexpr (isString) {
            if (*p_args == ArgType::s) {
                size_t strLength = std::strlen(head);
                buffer.write(strLength);
                buffer.writeBytes(head, strLength);
            }
            else
                buffer.write(head);
        }
        else
            buffer.write(head);
        writeArgs(p_args + 1, args...);
    }
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
