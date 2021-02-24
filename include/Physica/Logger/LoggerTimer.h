/*
 * Copyright 2020-2021 WeiBo He.
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

#include <ctime>
#include <cstdint>
#include "Physica/Utils/Cycler.h"

namespace Physica::Logger {
    /**
     * This class handle conversion between cycle at initialize and Unix time.
     */
    class LoggerTimer {
        uint64_t startCycle;
        time_t startTime;
    public:
        LoggerTimer();
        ~LoggerTimer() = default;
        /* Operations */
        [[nodiscard]] time_t now() const { return toTime(Utils::Cycler::now()); }
        [[nodiscard]] time_t toTime(uint64_t cycle) const;
        /* Getters */
        [[nodiscard]] uint64_t getStartCycle() const noexcept { return startCycle; }
        [[nodiscard]] time_t getStartTime() const noexcept { return startTime; }
    };
}
