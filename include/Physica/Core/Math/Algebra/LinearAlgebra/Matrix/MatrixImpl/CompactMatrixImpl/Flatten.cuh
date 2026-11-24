/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../CompactMatrix.cuh"
#include "Flatten.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    class device_obj<Flatten<M>> : public device_obj<CompactVector<Flatten<M>>> {
        using host_obj = Flatten<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactVector<host_obj>>;

        device_obj<M>& mat;
    public:
        device_obj(device_obj<CompactMatrix<M>>& mat_) : mat(mat_.getDerived()) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&&) noexcept;
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ auto device_obj<Flatten<M>>::data_handle(this auto&& self) noexcept {
        return self.mat.data_handle();
    }
}
