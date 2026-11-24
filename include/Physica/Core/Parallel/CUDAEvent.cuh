/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Exception/CUDA/CUDA.cuh"

namespace Physica {
    class PHYSICA_API CUDAEvent {
        using This = CUDAEvent;

        cudaEvent_t event = nullptr;
    public:
        enum : char {
            BlockingSync = cudaEventBlockingSync,
            DisableTiming = cudaEventDisableTiming,
            Interprocess = cudaEventInterprocess,
        };
    public:
        explicit CUDAEvent(int flag = cudaEventDefault);
        CUDAEvent(const This&) = delete;
        CUDAEvent(This&& obj) noexcept;
        ~CUDAEvent();
        /* Operators */
        CUDAEvent& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] operator cudaEvent_t() const noexcept { return event; }
        /* Operations */
        void record();
        inline void wait();
        void swap(CUDAEvent& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] cudaEvent_t getEvent() const noexcept { return event; }
    };

    inline void CUDAEvent::wait() {
        check(cudaEventSynchronize(event));
    }
}
