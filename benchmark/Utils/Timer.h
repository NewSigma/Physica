/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
#else
    class DeviceTimer;
#endif
}
