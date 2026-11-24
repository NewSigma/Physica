/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Macro.h"

struct CUstream_st;
using cudaStream_t = CUstream_st*;

namespace Physica {
    class PHYSICA_API CUDAStream {
        using This = CUDAStream;

        cudaStream_t stream = nullptr;
    public:
        CUDAStream();
        CUDAStream(std::nullptr_t);
        CUDAStream(const This&) = delete;
        CUDAStream(This&& obj) noexcept;
        ~CUDAStream();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        [[nodiscard]] __host__ __device__ operator cudaStream_t() const noexcept;
        /* Operations */
        void wait() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] bool query() const noexcept;
    };
}
