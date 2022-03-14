/*
 * Copyright 2020-2022 WeiBo He.
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

#include <cstdint>
#include <x86intrin.h>
#include "Physica/Core/Exception/RdrandException.h"

namespace Physica::Utils {
    class Random {
    public:
        static int retryLimit;
        /**
         * rdrand returns a high quality random number, in rare conditions it may fail and we should retry
         * several times but in extremely rare conditions it can not return anything.
         */
        static void rdrand(uint16_t& integer) {
            for(int i = 0; i < retryLimit; ++i) {
                const int code = _rdrand16_step(&integer);
                if (code == 1)
                    return;
            }
            throw Core::RdrandException();
        }

        static void rdrand(uint32_t& integer) {
            for(int i = 0; i < retryLimit; ++i) {
                const int code = _rdrand32_step(&integer);
                if (code == 1)
                    return;
            }
            throw Core::RdrandException();
        }

        static void rdrand(uint64_t& integer) {
            for(int i = 0; i < retryLimit; ++i) {
                const int code = _rdrand64_step(reinterpret_cast<unsigned long long*>(&integer));
                if (code == 1)
                    return;
            }
            throw Core::RdrandException();
        }
    };
}
