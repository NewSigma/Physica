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
#include <cassert>
#include <sys/time.h>
#include "Physica/Logger/LoggerTimer.h"

using namespace Physica;

LoggerTimer::LoggerTimer() {
    gettimeofday(&startTime, nullptr);
    startCycle = Cycler::now();
}

timeval LoggerTimer::toTime(uint64_t cycle) const {
    assert(cycle > startCycle);
    const uint64_t delta = cycle - startCycle;
    const uint64_t us = Cycler::toMicroseconds(delta) + startTime.tv_usec;
    const uint64_t s = us / 1000000;
    timeval result = startTime;
    result.tv_sec += s;
    result.tv_usec = us - 1000000 * s;
    return result;
}
