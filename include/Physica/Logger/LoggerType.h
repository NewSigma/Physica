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
#ifndef PHYSICA_LOGLEVEL_H
#define PHYSICA_LOGLEVEL_H

#include <cstdint>
#include <cstddef>
#include <array>

namespace Physica::Logger {
    /*!
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
        Fatal,
        Warning,
        Info,
        Debug,
        Off,
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

#endif
