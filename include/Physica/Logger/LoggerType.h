/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <cstdint>
#include <cstddef>
#include <array>

namespace Physica {
    /**
     * ArgType indicates the conversion specifier, e.g. c in %c, s in %s.
     */
    enum ArgType {
        c,
        s,
        d,
        o,
        x,
        X,
        u,
        U,
        f,
        F,
        a,
        A,
        g,
        G,
        p,
        Invalid
    };

    enum class LogLevel {
        Off,
        Fatal,
        Warn,
        Info,
        Debug,
        Global //Use global level instead.
    };

    struct LogInfo {
        const char* __restrict level;
        const char* __restrict format;
        const char* __restrict file;
        unsigned int line;
        const ArgType* args;
        size_t argCount;
    };
}
