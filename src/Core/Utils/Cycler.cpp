/*
 * Copyright 2020-2024 Weibo He.
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
#ifdef __linux__
    #include <sys/time.h>
#endif
#include <cassert>
#include <Physica/Core/Utils/Cycler.h>
#include <Physica/Core/Exception/SystemException.h>

namespace Physica::Core {
    /*!
     * Given an elapsed time measured in cycles, return an integer
     * giving the corresponding time in nanoseconds. Note: toSeconds()
     * is faster than this method.
     * \param cycles
     *      Difference between the results of two calls to rdtsc.
     * \param cyclesPerSec
     *      Optional parameter to specify the frequency of the counter that #cycles
     *      was taken from. Useful when converting a remote machine's tick counter
     *      to seconds. The default value of 0 will use the local processor's
     *      computed counter frequency.
     * \return
     *      The time in nanoseconds corresponding to cycles (rounded).
     */
    uint64_t Cycler::toNanoseconds(uint64_t cycles, double cyclesPerSec) {
        assert(cyclesPerSec > 0);
        if (cyclesPerSec == 0)
            cyclesPerSec = getCyclesPerSec();
        /*
         * Add 0.5 in order that when round from double to integer,
         * round to nearest integer instead of 0.
         * Avoid using library function lround() for performance.
         */
        return (uint64_t)(1e09 * static_cast<double>(cycles) / cyclesPerSec + 0.5); //NOLINT cycles > 0, cyclesPerSec > 0, rounding is safe
    }
    /*!
     * Given an elapsed time measured in cycles, return an integer
     * giving the corresponding time in microseconds. Note: toSeconds()
     * is faster than this method.
     * \param cycles
     *      Difference between the results of two calls to rdtsc.
     * \param cyclesPerSec
     *      Optional parameter to specify the frequency of the counter that #cycles
     *      was taken from. Useful when converting a remote machine's tick counter
     *      to seconds. The default value of 0 will use the local processor's
     *      computed counter frequency.
     * \return
     *      The time in microseconds corresponding to cycles (rounded).
     */
    uint64_t Cycler::toMicroseconds(uint64_t cycles, double cyclesPerSec) {
        assert(cyclesPerSec > 0);
        return toNanoseconds(cycles, cyclesPerSec) / 1000;
    }

    uint64_t Cycler::toMillisecond(uint64_t cycles, double cyclesPerSec) {
        assert(cyclesPerSec > 0);
        return toNanoseconds(cycles, cyclesPerSec) / 1000000;
    }
    /*!
     * Given an elapsed time measured in cycles, return a floating-point number
     * giving the corresponding time in seconds.
     * \param cycles
     *      Difference between the results of two calls to rdtsc.
     * \param cyclesPerSec
     *      Optional parameter to specify the frequency of the counter that #cycles
     *      was taken from. Useful when converting a remote machine's tick counter
     *      to seconds. The default value of 0 will use the local processor's
     *      computed counter frequency.
     * \return
     *      The time in seconds corresponding to cycles.
     */
    double Cycler::toSeconds(uint64_t cycles, double cyclesPerSec) {
        assert(cyclesPerSec > 0);
        if (cyclesPerSec == 0)
            cyclesPerSec = getCyclesPerSec();
        return static_cast<double>(cycles) / cyclesPerSec;
    }
    /*!
     * Given a number of nanoseconds, return an approximate number of
     * cycles for an equivalent time length.
     * \param ns
     *      Number of nanoseconds.
     * \param cyclesPerSec
     *      Optional parameter to specify the frequency of the counter that #cycles
     *      was taken from. Useful when converting a remote machine's tick counter
     *      to seconds. The default value of 0 will use the local processor's
     *      computed counter frequency.
     * \return
     *      The approximate number of cycles for the same time length.
     */
    uint64_t Cycler::fromNanoseconds(uint64_t ns, double cyclesPerSec) {
        assert(cyclesPerSec > 0);
        if (cyclesPerSec == 0)
            cyclesPerSec = getCyclesPerSec();
        /*
         * Add 0.5 in order that when round from double to integer,
         * round to nearest integer instead of 0.
         * Avoid using library function lround() for performance.
         */
        return (uint64_t)(static_cast<double>(ns) * cyclesPerSec / 1e09 + 0.5); //NOLINT ns > 0, cyclesPerSec > 0, rounding is safe
    }
    /*!
     * Given a time in seconds, return the number of cycles that it
     * corresponds to.
     * \param seconds
     *      Time in seconds.
     * \param cyclesPerSec
     *      Optional parameter to specify the frequency of the counter that #cycles
     *      was taken from. Useful when converting a remote machine's tick counter
     *      to seconds. The default value of 0 will use the local processor's
     *      computed counter frequency.
     * \return
     *      The approximate number of cycles corresponding to #seconds.
     */
    uint64_t Cycler::fromSeconds(double seconds, double cyclesPerSec) {
        assert(cyclesPerSec > 0);
        if (cyclesPerSec == 0)
            cyclesPerSec = getCyclesPerSec();
        /*
         * Add 0.5 in order that when round from double to integer,
         * round to nearest integer instead of 0.
         * Avoid using library function lround() for performance.
         */
        return (uint64_t)(seconds * cyclesPerSec + 0.5); //NOLINT second > 0, cyclesPerSec > 0, rounding is safe
    }

    double Cycler::makeCyclesPerSec() {
    #ifdef __linux__
        struct timeval startTime, stopTime; /* NOLINT no initializing is intended */
        // There is one tricky aspect, which is that we could get interrupted
        // between calling gettimeofday and reading the cycle counter, in which
        // case we won't have corresponding readings.  To handle this (unlikely)
        // case, compute the overall result repeatedly, and wait until we get
        // two successive calculations that are within 0.001% of each other (or
        // in other words, a drift of up to 10µs per second).
        double oldCycles = 0;
        double result = 0;
        while (true) {
            oldCycles = result;
            if (gettimeofday(&startTime, nullptr) != 0) [[unlikely]]
                throw Core::SystemException();

            const uint64_t startCycles = now();
            while (true) {
                if (gettimeofday(&stopTime, nullptr) != 0) [[unlikely]]
                    throw Core::SystemException();

                const uint64_t stopCycles = now();
                const uint64_t micros = (stopTime.tv_usec - startTime.tv_usec) +
                                        (stopTime.tv_sec - startTime.tv_sec) * 1000000;
                if (micros > 10000) {
                    result = static_cast<double>(stopCycles - startCycles);
                    result = 1000000.0 * result / static_cast<double>(micros);
                    break;
                }
            }
            const double delta = result / 100000.0;
            if (((result - delta) < oldCycles) && (oldCycles < (result + delta)))
                break;
        }
        return result;
    #else
        noImpl();
    #endif
    }
}
