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
#ifndef PHYSICA_CYCLER_H
#define PHYSICA_CYCLER_H

#include <cstdint>

namespace Physica::Utils {
    /**
     * This class provides static methods that read the fine-grain CPU
     * cycle counter and translate between cycle-level times and absolute
     * times.
     *
     * Warning:
     * It is rare that some operating system will not enable the users use rdstc instruction,
     * but they are likely to provide a interface for users to call.
     *
     * Questions:
     * 1.According to [1], we can force a program run in kernel to avoid interrupts, may be we can use of it?
     * 2.According to [1], we can make use of asm CPUID to get a preciser result.
     *
     * Reference:
     * [1] PAOLONI, G. How to benchmark code execution times on intel ia-32 and ia-64 instruction set
     * architectures. Intel Corporation, September 123 (2010).
     * [2] NanoLog (https://github.com/PlatformLab/NanoLog)
     */
    class Cycler {
        /*!
         * Conversion factor between cycles and the seconds; computed by Cycler::init.
         */
        static double localCyclesPerSec;
    public:
        /* Static Members */
        static void init();
        /**
         * Return the current value of the fine-grain CPU cycle counter
         *
         * Because of out of order execution, the result may be less  or more than several cycles.
         * Use tic() and toc() to avoid this problem.
         */
        static __inline __attribute__((always_inline))
        uint64_t now() {
            uint32_t lo, hi;
            __asm__ __volatile__ (
                    "rdtsc"
                    : "=a" (lo), "=d" (hi)
                    ::
            );
            return ((static_cast<uint64_t>(hi) << 32U) | lo);
        }
        /**
         * Start the timing, the instructions
         * before this function is all executed, Use this function and toc() in pair.
         *
         * This function is slower than now().
         */
        static __inline __attribute__((always_inline))
        uint64_t tic() {
            uint32_t lo, hi;
            __asm__ __volatile__ (
                    "cpuid\n\t"
                    "rdtsc\n\t"
                    : "=a" (lo), "=d" (hi)
                    :: "ebx", "ecx"
            );
            return ((static_cast<uint64_t>(hi) << 32U) | lo);
        }
        /**
         * End the timing, the instructions
         * before this function is all executed, Use this function and tic() in pair.
         *
         * This function is slower than now().
         */
        static __inline __attribute__((always_inline))
        uint64_t toc() {
            uint32_t lo, hi;
            __asm__ __volatile__ (
                    "rdtscp\n\t"
                    : "=a" (lo), "=d" (hi)
            );
            __asm__ __volatile__ (
                    "cpuid\n\t"
                    ::: "eax", "ebx", "ecx", "edx"
            );
            return ((static_cast<uint64_t>(hi) << 32U) | lo);
        }
        /**
         * Returns the conversion factor between cycles in seconds, using
         * a mock value for testing when appropriate.
         */
        [[nodiscard]] static __inline __attribute__((always_inline))
        double getCyclesPerSec() { return localCyclesPerSec; }

        static uint64_t toNanoseconds(uint64_t cycles, double cyclesPerSec = localCyclesPerSec);
        static uint64_t toMicroseconds(uint64_t cycles, double cyclesPerSec = localCyclesPerSec);
        static double toSeconds(uint64_t cycles, double cyclesPerSec = localCyclesPerSec);
        static uint64_t fromNanoseconds(uint64_t ns, double cyclesPerSec = localCyclesPerSec);
        static uint64_t fromSeconds(double seconds, double cyclesPerSec = localCyclesPerSec);
    private:
        Cycler() = default;
    };
}

#endif
