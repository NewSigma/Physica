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
#ifndef PHYSICA_FORMATANALYZERIMPL_H
#define PHYSICA_FORMATANALYZERIMPL_H

#include <utility>

namespace Physica::Logger {
    /*!
     * \param format
     *      C style format string.
     *
     * \return
     *      Return the dynamic arguments in format.
     */
    template<size_t N>
    constexpr size_t getArgCount(const char (&format)[N]) {
        size_t pos = 0;
        size_t argCount = 0;
        while(pos < N - 1) {
            if(format[pos] == '%' && format[pos + 1] != '%')
                ++argCount;
            ++pos;
        }
        return argCount;
    }
    /*!
     * Get the arg type at index position.
     *
     * \tparam N
     *      Length of the static format string (automatically deduced)
     * \param format
     *      Format string to parse
     * \param index
     *      The index of arg type, starts from 0.
     * \return
     *      Returns the arg type at index position.
     */
    template<size_t N>
    constexpr ArgType getArgType(const char (&format)[N], int index) {
        size_t pos = 0;
        while (pos < N - 1) {
            if (format[pos] != '%') {
                ++pos;
                continue;
            }
            else {
                // Note: gcc++ 5,6,7,8 seems to hang whenever one uses the construct
                // "if (...) {... continue; }" without an else in constexpr
                // functions. Hence, we have the code here wrapped in an else {...}
                // Reference:
                // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86767
                ++pos;

                // Two %'s in a row => Comment
                if (format[pos] == '%') {
                    ++pos;
                    continue;
                }
                else {
                    ArgType result = ArgType::Invalid;
                    switch (format[pos]) {
                        case 'c':
                            result = c;
                            break;
                        case 's':
                            result = s;
                            break;
                        case 'd':
                        case 'i':
                            result = d;
                            break;
                        case 'o':
                            result = o;
                            break;
                        case 'x':
                            result = x;
                            break;
                        case 'X':
                            result = X;
                            break;
                        case 'u':
                            result = u;
                            break;
                        case 'f':
                            result = f;
                            break;
                        case 'F':
                            result = F;
                            break;
                        case 'a':
                            result = a;
                            break;
                        case 'A':
                            result = A;
                            break;
                        case 'g':
                            result = g;
                            break;
                        case 'G':
                            result = G;
                            break;
                        case 'p':
                            result = p;
                            break;
                        case 'n':
                            throw std::invalid_argument("%n specifier are not support!");
                        default:
                            throw std::invalid_argument("Unrecognized format specifier after %");
                    }
                    if(index == 0)
                        return result;
                    else
                        --index;
                }
            }
        }
        return ArgType::Invalid;
    }
    /*!
     * Helper to analyzeFormatString. This level of indirection is needed to
     * unpack the index_sequence generated in analyzeFormatString and
     * use the sequence as indices for calling getArgType.
     *
     * \tparam N
     *      Length of the format string (automatically deduced)
     * \tparam Indices
     *      An index sequence from [0, N) where N is the number of parameters in
     *      the format string (automatically deduced)
     * \param format
     *      printf format string to analyze
     * \return
     *      An std::array describing the types at each index (zero based).
     */
    template<size_t N, std::size_t... Indices>
    constexpr std::array<ArgType, sizeof...(Indices)>
    analyzeFormatStringHelper(const char(&format)[N], std::index_sequence<Indices...>) {
        return {{ getArgType(format, Indices)... }};
    }
    /*!
     * Spawn a array of all of the args of a format string.
     *
     * \tparam ArgCount
     *      The number of additional format parameters that follow the format
     *      string in a printf-like function. For example printf("%*.*d", 9, 8, 7)
     *      would have NParams = 3
     * \tparam N
     *      length of the printf style format string (automatically deduced)
     * \param format
     *      Format string to generate the array for
     * \return
     *      An std::array where the n-th index indicates a conversion specifier.
     */
    template<size_t ArgCount, size_t N>
    constexpr std::array<ArgType, ArgCount> analyzeFormatString(const char(&format)[N]) {
        return analyzeFormatStringHelper(format, std::make_index_sequence<ArgCount>{});
    }
}

#endif
