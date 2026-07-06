/*
 * Copyright 2023 Weibo He.
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

namespace Physica {
    /**
     * Reference:
     * [1] vectorclass2; https://github.com/vectorclass/version2/blob/master/instrset.h
     */
    class Instruset {
    public:
        [[nodiscard]] consteval static bool hasAVX512VL() {
        #ifdef __AVX512VL__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasAVX512BW() {
        #ifdef __AVX512BW__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasAVX512DQ() {
        #ifdef __AVX512DQ__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasAVX512() {
        #if defined (__AVX512F__) || defined (__AVX512__)
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasAVX2() {
            if constexpr (hasAVX512())
                return true;
        #ifdef __AVX2__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasAVX() {
            if constexpr (hasAVX2())
                return true;
        #ifdef __AVX__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasSSE4_2() {
            if constexpr (hasAVX())
                return true;
        #ifdef __SSE_4_2__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasSSE4_1() {
            if constexpr (hasSSE4_2())
                return true;
        #ifdef __SSE_4_1__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasSSSE3() {
            if constexpr (hasSSE4_1())
                return true;
        #ifdef __SSSE3__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasSSE3() {
            if constexpr (hasSSSE3())
                return true;
        #ifdef __SSE3__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasSSE2() {
            if constexpr (hasSSE3())
                return true;
        #if defined (__SSE2__) || defined (__x86_64__)
            return true;
        #elif defined _M_IX86_FP
            return _M_IX86_FP == 2;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasSSE() {
            if constexpr (hasSSE2())
                return true;
        #ifdef __SSE__
            return true;
        #elif defined _M_IX86_FP
            return _M_IX86_FP == 1;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasFMA() {
        #ifdef __FMA__
            return true;
        #else
            return false;
        #endif
        }

        [[nodiscard]] consteval static bool hasFMA4() {
        #ifdef __FMA4__
            return true;
        #else
            return false;
        #endif
        }
    };
}
