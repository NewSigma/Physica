/*
 * Copyright 2020-2026 Weibo He.
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
#include <cmath>
#include <iostream>
#ifdef __linux__
    #include <sys/time.h>
#endif
#include <cassert>
#include "Physica/Core/Utils/Cycler.h"

using namespace Physica;

namespace {
    timeval getTimeOfDay() noexcept {
        timeval result{};
        if (gettimeofday(&result, nullptr) != 0) [[unlikely]] {
            std::cerr << "Cycler: gettimeofday failed";
            std::abort();
        }
        return result;
    }

    double makeCyclesPerSec() noexcept {
    #ifdef __linux__
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
            const auto startTime = getTimeOfDay();

            const uint64_t startCycles = Cycler::now();
            while (true) {
                const auto stopTime = getTimeOfDay();
                const uint64_t stopCycles = Cycler::now();
                const uint64_t micros = (stopTime.tv_usec - startTime.tv_usec) + (stopTime.tv_sec - startTime.tv_sec) * 1000000;
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

double Cycler::getCyclesPerSec() noexcept {
    static const double localCyclesPerSec = makeCyclesPerSec();
    return localCyclesPerSec;
}

uint64_t Cycler::toNanoseconds(uint64_t cycles, double cyclesPerSec) noexcept {
    assert(cyclesPerSec > 0);
    if (cyclesPerSec == 0)
        cyclesPerSec = getCyclesPerSec();
    return std::llround(1e09 * static_cast<double>(cycles) / cyclesPerSec);
}

uint64_t Cycler::toMicroseconds(uint64_t cycles, double cyclesPerSec) noexcept {
    assert(cyclesPerSec > 0);
    return toNanoseconds(cycles, cyclesPerSec) / 1000;
}

uint64_t Cycler::toMillisecond(uint64_t cycles, double cyclesPerSec) noexcept {
    assert(cyclesPerSec > 0);
    return toNanoseconds(cycles, cyclesPerSec) / 1000000;
}

double Cycler::toSeconds(uint64_t cycles, double cyclesPerSec) noexcept {
    assert(cyclesPerSec > 0);
    if (cyclesPerSec == 0)
        cyclesPerSec = getCyclesPerSec();
    return static_cast<double>(cycles) / cyclesPerSec;
}

uint64_t Cycler::fromNanoseconds(uint64_t ns, double cyclesPerSec) noexcept {
    assert(cyclesPerSec > 0);
    if (cyclesPerSec == 0)
        cyclesPerSec = getCyclesPerSec();
    return std::llround(static_cast<double>(ns) * cyclesPerSec / 1e09);
}

uint64_t Cycler::fromSeconds(double seconds, double cyclesPerSec) noexcept {
    assert(cyclesPerSec > 0);
    if (cyclesPerSec == 0)
        cyclesPerSec = getCyclesPerSec();
    return std::llround(seconds * cyclesPerSec);
}
