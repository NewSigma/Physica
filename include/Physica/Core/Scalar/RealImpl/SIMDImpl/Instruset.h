/*
 * Copyright 2023-2026 Weibo He.
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
    class Instruset final {
    public:
        [[nodiscard]] consteval static bool hasAVX512VL() noexcept;
        [[nodiscard]] consteval static bool hasAVX512BW() noexcept;
        [[nodiscard]] consteval static bool hasAVX512DQ() noexcept;
        [[nodiscard]] consteval static bool hasAVX512() noexcept;
        [[nodiscard]] consteval static bool hasAVX2() noexcept;
        [[nodiscard]] consteval static bool hasAVX() noexcept;
        [[nodiscard]] consteval static bool hasSSE4_2() noexcept;
        [[nodiscard]] consteval static bool hasSSE4_1() noexcept;
        [[nodiscard]] consteval static bool hasSSSE3() noexcept;
        [[nodiscard]] consteval static bool hasSSE3() noexcept;
        [[nodiscard]] consteval static bool hasSSE2() noexcept;
        [[nodiscard]] consteval static bool hasSSE() noexcept;
        [[nodiscard]] consteval static bool hasFMA() noexcept;
        [[nodiscard]] consteval static bool hasFMA4() noexcept;
        /* Aliases */
        [[nodiscard]] consteval static bool support128() noexcept { return hasSSE2(); }
        [[nodiscard]] consteval static bool support256() noexcept { return hasSSE4_2(); }
        [[nodiscard]] consteval static bool support512() noexcept { return hasAVX512(); }
    };

    consteval bool Instruset::hasAVX512VL() noexcept {
    #ifdef __AVX512VL__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasAVX512BW() noexcept {
    #ifdef __AVX512BW__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasAVX512DQ() noexcept {
    #ifdef __AVX512DQ__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasAVX512() noexcept {
    #if defined (__AVX512F__) || defined (__AVX512__)
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasAVX2() noexcept {
        if constexpr (hasAVX512())
            return true;
    #ifdef __AVX2__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasAVX() noexcept {
        if constexpr (hasAVX2())
            return true;
    #ifdef __AVX__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasSSE4_2() noexcept {
        if constexpr (hasAVX())
            return true;
    #ifdef __SSE_4_2__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasSSE4_1() noexcept {
        if constexpr (hasSSE4_2())
            return true;
    #ifdef __SSE_4_1__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasSSSE3() noexcept {
        if constexpr (hasSSE4_1())
            return true;
    #ifdef __SSSE3__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasSSE3() noexcept {
        if constexpr (hasSSSE3())
            return true;
    #ifdef __SSE3__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasSSE2() noexcept {
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

    consteval bool Instruset::hasSSE() noexcept {
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

    consteval bool Instruset::hasFMA() noexcept {
    #ifdef __FMA__
        return true;
    #else
        return false;
    #endif
    }

    consteval bool Instruset::hasFMA4() noexcept {
    #ifdef __FMA4__
        return true;
    #else
        return false;
    #endif
    }
}
