/*
 * Copyright 2020-2025 Weibo He.
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

#ifdef __linux__
    #include <x86intrin.h>
#endif
#include <cstdint>
#include "Physica/Core/Exception/RdrandException.h"

namespace Physica {
    /**
     * rdrand returns a high quality random number, in rare conditions it may fail and we should retry
     * several times but in extremely rare conditions it can not return anything.
     */
    class RandomSeed {
        constexpr static int RetryLimit = 8;
    public:
    #ifdef __linux__
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
        /**
         * Generate a sequence of random seed (using the PCG-XSH-RS scheme)
         *
         * Reference:
         * [1] Eigen; https://eigen.tuxfamily.org/  
         * [2] https://www.pcg-random.org/
         */
        static uint32_t toNextSeed(uint64_t& state) noexcept {
            const uint64_t current = state;
            state = current * 6364136223846793005ULL + 0xda3e39cb94b95bdbULL;
            return static_cast<uint32_t>((current ^ (current >> 22U)) >> (22U + (current >> 61U)));
        }
    };
}
