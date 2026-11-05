/*
 * Copyright 2025-2026 Weibo He.
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
#include <tuple>
#include "Physica/Core/Utils/Cycler.h"
#ifdef PHYSICA_CUDA
    #include "Physica/Core/Parallel/CUDAEvent.cuh"
#endif

namespace Physica {
    class HostTimer {
        using This = HostTimer;

        uint64_t from{};
    public:
        HostTimer() noexcept;
        HostTimer(const This&) = delete;
        HostTimer(This&&) noexcept = delete;
        ~HostTimer() = default;
        /* Operators */
        HostTimer& operator=(const This&) = delete;
        HostTimer& operator=(This&&) noexcept = delete;
        /* Operations */
        [[gnu::always_inline]] void tic() noexcept;
        [[nodiscard, gnu::always_inline]] double toc() const noexcept;
    };

    inline HostTimer::HostTimer() noexcept {
        std::ignore = Cycler::getCyclesPerSec();
    }

    [[gnu::always_inline]] inline void HostTimer::tic() noexcept {
        from = Cycler::tic();
    }

    [[nodiscard, gnu::always_inline]] inline double HostTimer::toc() const noexcept {
        return Cycler::toSeconds(Cycler::toc() - from);
    }

#ifdef PHYSICA_CUDA
    class DeviceTimer {
        using This = DeviceTimer;

        CUDAEvent from;
        CUDAEvent to;
    public:
        DeviceTimer() = default;
        DeviceTimer(const This&) = delete;
        DeviceTimer(This&&) noexcept = delete;
        ~DeviceTimer() = default;
        /* Operators */
        DeviceTimer& operator=(const This&) = delete;
        DeviceTimer& operator=(This&&) noexcept = delete;
        /* Operations */
        [[gnu::always_inline]] void tic();
        [[nodiscard, gnu::always_inline]] double toc();
    };

    [[gnu::always_inline]] inline void DeviceTimer::tic() {
        from.record();
    }

    [[nodiscard, gnu::always_inline]] inline double DeviceTimer::toc() {
        to.record();
        to.wait();
        float elapsedMs = 0.0F;
        check(cudaEventElapsedTime(&elapsedMs, from, to));
        return static_cast<double>(elapsedMs) * 1E-3;
    }
#endif
}
