/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <ctime>
#include <cstdint>
#include "Physica/Core/Utils/Cycler.h"

namespace Physica {
    /**
     * This class handle conversion between cycle at initialize and Unix time.
     */
    class PHYSICA_API LoggerTimer {
        timeval startTime;
        uint64_t startCycle;
    public:
        LoggerTimer();
        ~LoggerTimer() = default;
        /* Operations */
        [[nodiscard]] timeval now() const { return toTime(Cycler::now()); }
        [[nodiscard]] timeval toTime(uint64_t cycle) const;
        /* Getters */
        [[nodiscard]] timeval getStartTime() const noexcept { return startTime; }
        [[nodiscard]] uint64_t getStartCycle() const noexcept { return startCycle; }
    };
}
