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
    /*!
     * Describes the type of parameter that would be passed into a printf-like
     * function.
     *
     * These types are optimized to store enough information to determine
     * (a) whether a 'const char*' parameter indicates string (%s) or not (%p)
     * (b) if a string parameter (%s) needs to be truncated due to precision
     * (c) whether a parameter is a dynamic precision/width specifier
     */
    enum ParamType : int32_t {
        // Indicates that there is a problem with the parameter
        INVALID = -6,

        // Indicates a dynamic width (i.e. the '*' in  %*.d)
        DYNAMIC_WIDTH = -5,

        // Indicates dynamic precision (i.e. the '*' in %.*d)
        DYNAMIC_PRECISION = -4,

        // Indicates that the parameter is not a string type (i.e. %d, %lf)
        NON_STRING = -3,

        // Indicates the parameter is a string and has a dynamic precision
        // (i.e. '%.*s' )
        STRING_WITH_DYNAMIC_PRECISION = -2,

        // Indicates a string with no precision specified (i.e. '%s' )
        STRING_WITH_NO_PRECISION = -1,

        // All non-negative values indicate a string with a precision equal to its
        // enum value casted as an int32_t
        STRING = 0
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
        int line;
    };
}

#endif
