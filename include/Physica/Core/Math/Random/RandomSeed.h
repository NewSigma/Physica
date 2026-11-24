/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <exception>
#include <cstdint>
#ifdef __linux__
    #include <x86intrin.h>
#endif

namespace Physica {
    /**
     * rdrand returns a high quality random number, in rare conditions it may fail and we should retry
     * several times but in extremely rare conditions it can not return anything.
     */
    class RandomSeed {
        constexpr static int RetryLimit = 8;

        class RdrandException : public std::exception {};
    public:
    #ifdef __linux__
        static void rdseed(uint16_t& __restrict integer) {
            for(int i = 0; i < RetryLimit; ++i) {
                const int code = _rdseed16_step(&integer);
                if (code == 1) [[likely]]
                    return;
            }
            throw RdrandException();
        }

        static void rdseed(uint32_t& __restrict integer) {
            for(int i = 0; i < RetryLimit; ++i) {
                const int code = _rdseed32_step(&integer);
                if (code == 1) [[likely]]
                    return;
            }
            throw RdrandException();
        }

        static void rdseed(uint64_t& __restrict integer) {
            for(int i = 0; i < RetryLimit; ++i) {
                const int code = _rdseed64_step(reinterpret_cast<unsigned long long*>(&integer));
                if (code == 1) [[likely]]
                    return;
            }
            throw RdrandException();
        }

        static void rdrand(uint16_t& __restrict integer) {
            for(int i = 0; i < RetryLimit; ++i) {
                const int code = _rdrand16_step(&integer);
                if (code == 1) [[likely]]
                    return;
            }
            throw RdrandException();
        }

        static void rdrand(uint32_t& __restrict integer) {
            for(int i = 0; i < RetryLimit; ++i) {
                const int code = _rdrand32_step(&integer);
                if (code == 1) [[likely]]
                    return;
            }
            throw RdrandException();
        }

        static void rdrand(uint64_t& __restrict integer) {
            for(int i = 0; i < RetryLimit; ++i) {
                const int code = _rdrand64_step(reinterpret_cast<unsigned long long*>(&integer));
                if (code == 1) [[likely]]
                    return;
            }
            throw RdrandException();
        }
    #endif
    };
}
