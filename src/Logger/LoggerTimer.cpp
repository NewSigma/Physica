/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
