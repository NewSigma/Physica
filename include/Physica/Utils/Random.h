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
#ifndef PHYSICA_RANDOM_H
#define PHYSICA_RANDOM_H

#include <x86intrin.h>

namespace Physica::Utils {
    class Random {
        static unsigned char retryLimit;
        /*!
         * rdrand returns a high quality random number, in rare conditions it may fail and we should retry
         * several times but in extremely rare conditions it can not return anything.
         */
        static void rdrand(__uint16_t& integer) {
            for(int i = 0; i < retryLimit; i ++) {
                if(_rdrand16_step(&integer))
                    break;
            }
            abort();
        }

        static void rdrand(__uint32_t& integer) {
            for(int i = 0; i < retryLimit; i ++) {
                if(_rdrand32_step(&integer))
                    break;
            }
            abort();
        }

        static void rdrand(__uint64_t& integer) {
            for(int i = 0; i < retryLimit; i ++) {
                if(_rdrand64_step(reinterpret_cast<unsigned long long*>(&integer)))
                    break;
            }
            abort();
        }
    };
}

#endif
